#!/bin/bash
# Prepare missing lanthanide directories, then submit a per-metal SLURM
# afterok chain. Existing directories are SKIPPED unless the user explicitly
# opts in to submitting prepared ones (--include-existing).
set -Eeuo pipefail

usage() {
cat <<'HELP'
Submit TURBOMOLE lanthanide jobs on VACC (14 metals; Pm excluded).
Run from an EF project directory containing control, coord and def.inp.

Usage: lan-run [num|ao] [options]       (default: no Hessian)
  lan-run = geometry optimization only (no Hessian)
  num = NumForce; ao = aoforce (full names also accepted).

Options:
  --metals sm tb      Only these metals (default: all 14)
  --prepare-only      Prepare folders; submit nothing
  --include-existing  Submit already-prepared folders not previously chain-submitted
  --base DIR          Master project directory (default: current directory)
  -h, --help          Show help

Examples:
  lan-run num --metals sm --prepare-only  # prepare Sm for inspection
  lan-run num --metals sm --include-existing  # submit prepared Sm
  lan-run                                 # create/submit missing metals, no Hessian
  lan-run ao                              # create/submit missing metals, aoforce
  lan-run num --metals gd tb              # two metals, NumForce

Existing metal folders are never modified. New gd/ requires master def_gd.inp.
Jobs are connected by afterok; orbital holes are warnings, not failures.
Help for difficult SCF cases: lan-recover --help
HELP
}

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
base="$PWD"
frequency='none'
prepare_only=0
include_existing=0
metals=(la ce pr nd sm eu gd tb dy ho er tm yb lu)
selected_metals=()

while (($#)); do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        num|numforce) frequency='numforce'; shift ;;
        ao|aoforce) frequency='aoforce'; shift ;;
        --base)
            (($# >= 2)) || { echo 'ERROR: --base needs a directory' >&2; exit 2; }
            base="$2"; shift 2 ;;
        --metals)
            shift
            selected_metals=()
            while (($#)) && [[ "$1" != -* && "$1" != numforce && "$1" != aoforce && "$1" != num && "$1" != ao ]]; do
                selected_metals+=("${1,,}")
                shift
            done
            ((${#selected_metals[@]})) || { echo 'ERROR: --metals requires at least one symbol' >&2; exit 2; }
            ;;
        --prepare-only) prepare_only=1; shift ;;
        --include-existing) include_existing=1; shift ;;
        *) echo "ERROR: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

if ((${#selected_metals[@]})); then metals=("${selected_metals[@]}"); fi
for metal in "${metals[@]}"; do
    case "$metal" in
        la|ce|pr|nd|sm|eu|gd|tb|dy|ho|er|tm|yb|lu) ;;
        *) echo "ERROR: invalid lanthanide '$metal' (Pm is intentionally excluded)" >&2; exit 2 ;;
    esac
done

base="$(cd -- "$base" && pwd -P)"
for input in control coord def.inp; do
    [[ -f "$base/$input" ]] || { echo "ERROR: master $input missing in $base" >&2; exit 1; }
done
for helper in setup_lanthanides.py tm_workflow.py turbo.pipeline.optimized.slurm; do
    [[ -f "$SCRIPT_DIR/$helper" ]] || { echo "ERROR: $helper missing from package" >&2; exit 1; }
done

: "${TM_CPUS:=20}" "${EIGER_CPUS:=1}" "${RIDFT_MEM:=20G}" "${JOBEX_MEM:=20G}" "${HESS_MEM:=100G}" "${EIGER_MEM:=4G}"
[[ "$TM_CPUS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid TM_CPUS=$TM_CPUS" >&2; exit 2; }
[[ "$EIGER_CPUS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid EIGER_CPUS=$EIGER_CPUS" >&2; exit 2; }

created=()
for metal in "${metals[@]}"; do
    if [[ -e "$base/$metal" ]]; then
        echo "SKIP SETUP $metal: existing path preserved completely"
        continue
    fi
    python3 "$SCRIPT_DIR/setup_lanthanides.py" --base "$base" --metals "$metal"
    created+=("$metal")
done

if ((prepare_only)); then
    echo "Preparation finished. Created ${#created[@]} folder(s); no jobs submitted."
    exit 0
fi

# Submit only directories this invocation created, unless deliberately opted in
# to PREPARED existing directories; never modify existing inputs.
submit_metals=("${created[@]}")
if ((include_existing)); then
    for metal in "${metals[@]}"; do
        for newmetal in "${created[@]}"; do
            [[ "$metal" != "$newmetal" ]] || continue 2
        done
        submit_metals+=("$metal")
    done
fi

submit_metal() {
    local metal="$1" dir="$base/$1" stage jid dep='' memory time_limit raw cpus partition
    if [[ ! -f "$dir/metal.json" || ! -f "$dir/control" || ! -f "$dir/coord" || \
          ! -f "$dir/def.inp" || ! -f "$dir/tm_workflow.py" || \
          ! -f "$dir/turbo.pipeline.optimized.slurm" || -e "$dir/control.template" ]]; then
        echo "SKIP SUBMIT $metal: not a valid prepared folder or control.template present"
        return 0
    fi
    if [[ -e "$dir/chain_submission.tsv" ]]; then
        echo "SKIP SUBMIT $metal: chain_submission.tsv exists (prevents duplicate jobs)"
        return 0
    fi

    # The submit call must run FROM the metal directory, so SLURM_SUBMIT_DIR
    # and all TURBOMOLE output files refer to that directory.
    local stages=(ridft eiger_pre occupation_fix integer_ridft jobex eiger_post)
    if [[ "$frequency" != none ]]; then stages+=("$frequency"); fi
    for stage in "${stages[@]}"; do
        partition='general'
        case "$stage" in
            eiger_*) time_limit='00:10:00'; memory="$EIGER_MEM"; cpus="$EIGER_CPUS" ;;
            occupation_fix) partition='short'; time_limit='00:10:00'; memory='1G'; cpus='1' ;;
            ridft|integer_ridft) time_limit='2-00:00:00'; memory="$RIDFT_MEM"; cpus="$TM_CPUS" ;;
            jobex) time_limit='2-00:00:00'; memory="$JOBEX_MEM"; cpus="$TM_CPUS" ;;
            numforce|aoforce) time_limit='2-00:00:00'; memory="$HESS_MEM"; cpus="$TM_CPUS" ;;
        esac
        local export_spec="ALL,TM_STAGE=$stage"
        if [[ "$frequency" == none && "$stage" == eiger_post ]]; then
            export_spec+=",TM_FINAL_NO_HESSIAN=1"
        fi
        local opts=(--parsable --partition="$partition" --time="$time_limit"
                    --mem="$memory" --cpus-per-task="$cpus"
                    --job-name="tm-${metal}-${stage}"
                    --export="$export_spec")
        if [[ -n "$dep" ]]; then opts+=("--dependency=afterok:$dep"); fi
        if ! raw="$(cd "$dir" && sbatch "${opts[@]}" turbo.pipeline.optimized.slurm)"; then
            echo "ERROR $metal: sbatch failed at $stage. Previously submitted stages are listed in $dir/chain_submission.tsv." >&2
            return 1
        fi
        jid="${raw%%;*}"
        [[ "$jid" =~ ^[0-9]+$ ]] || {
            echo "ERROR $metal: unexpected sbatch output for $stage: $raw" >&2
            return 1
        }
        if [[ -z "$dep" ]]; then
            printf 'stage\tjob_id\tpartition\twalltime\tmemory\n' > "$dir/chain_submission.tsv"
        fi
        printf '%s\t%s\t%s\t%s\t%s\n' "$stage" "$jid" "$partition" "$time_limit" "$memory" >> "$dir/chain_submission.tsv"
        printf '%-2s %-14s %-10s %s\n' "$metal" "$stage" "$jid" "${dep:+afterok:$dep}"
        dep="$jid"
    done
    if [[ "$frequency" == none ]]; then
        echo "SUBMITTED $metal: 6-stage geometry-only chain (no Hessian; final SLURM job $dep)"
    else
        echo "SUBMITTED $metal: 7-stage $frequency chain (final SLURM job $dep)"
    fi
}

for metal in "${submit_metals[@]}"; do submit_metal "$metal"; done
printf 'Done. Prepared %d new folder(s). See each metal/chain_submission.tsv.\n' "${#created[@]}"
echo 'To monitor: squeue -u "$USER"'
echo 'On failure, afterok blocks all dependent stages; orbital holes are warnings only.'
