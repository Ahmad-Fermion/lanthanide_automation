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
  --include-existing  Submit pristine, never-run folders from lan-setup
  --resume-partial    Resume an interrupted chain without repeating submitted jobs
                      (requires --metals and matching original Hessian selection)
  --base DIR          Master project directory (default: current directory)
  -h, --help          Show help

Examples:
  lan-run num --metals sm --prepare-only  # prepare Sm for inspection
  lan-run num --metals sm --include-existing  # submit pristine, prepared Sm
  lan-run num --metals sm --resume-partial    # resume missing submissions only
  lan-run                                 # create/submit missing metals, no Hessian
  lan-run ao                              # create/submit missing metals, aoforce
  lan-run num --metals gd tb              # two metals, NumForce

Existing scientific inputs are preserved; opt-in submissions add only job records.
New gd/ requires master def_gd.inp.
Jobs are connected by afterok; orbital holes are warnings, not failures.
For a partial submission, --resume-partial verifies the original chain plan,
SLURM job status and work directory. It never reruns DEFINE or RIDFT.
If a scientific stage failed, use lan-recover after inspecting the calculation.
Help for difficult SCF cases: lan-recover --help
HELP
}

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
base="$PWD"
frequency='none'
prepare_only=0
include_existing=0
resume_partial=0
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
        --resume-partial) resume_partial=1; shift ;;
        *) echo "ERROR: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

if ((resume_partial)); then
    (( ! include_existing && ! prepare_only )) || {
        echo 'ERROR: --resume-partial cannot be combined with --include-existing or --prepare-only' >&2; exit 2;
    }
    ((${#selected_metals[@]})) || {
        echo 'ERROR: --resume-partial requires explicit --metals (only known partial chains)' >&2; exit 2;
    }
fi
if ((${#selected_metals[@]})); then metals=("${selected_metals[@]}"); fi
for metal in "${metals[@]}"; do
    case "$metal" in
        la|ce|pr|nd|sm|eu|gd|tb|dy|ho|er|tm|yb|lu) ;;
        *) echo "ERROR: invalid lanthanide '$metal' (Pm is intentionally excluded)" >&2; exit 2 ;;
    esac
done

base="$(cd -- "$base" && pwd -P)"
if (( ! resume_partial )); then
    for input in control coord def.inp; do
        [[ -f "$base/$input" ]] || { echo "ERROR: master $input missing in $base" >&2; exit 1; }
    done
fi
for helper in setup_lanthanides.py submission_guard.py tm_workflow.py turbo.pipeline.optimized.slurm; do
    [[ -f "$SCRIPT_DIR/$helper" ]] || { echo "ERROR: $helper missing from package" >&2; exit 1; }
done

: "${TM_CPUS:=20}" "${EIGER_CPUS:=1}" "${RIDFT_MEM:=20G}" "${JOBEX_MEM:=20G}" "${HESS_MEM:=100G}" "${EIGER_MEM:=4G}"
[[ "$TM_CPUS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid TM_CPUS=$TM_CPUS" >&2; exit 2; }
[[ "$EIGER_CPUS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid EIGER_CPUS=$EIGER_CPUS" >&2; exit 2; }

created=()
failures=0
if (( ! resume_partial )); then
    for metal in "${metals[@]}"; do
        if [[ -e "$base/$metal" ]]; then
            echo "SKIP SETUP $metal: existing path preserved completely"
            continue
        fi
        if python3 "$SCRIPT_DIR/setup_lanthanides.py" --base "$base" --metals "$metal"; then
            created+=("$metal")
        else
            echo "ERROR $metal: preparation failed; continuing with other metals" >&2
            ((failures+=1))
        fi
    done
fi

if ((prepare_only)); then
    echo "Preparation finished. Created ${#created[@]} folder(s); no jobs submitted. Failed: $failures."
    ((failures == 0))
    exit $?
fi

# Normal operation submits newly created folders and (if explicitly requested)
# pristine, previously prepared folders. Partial-chain continuation is NEVER
# automatic: the user must explicitly request --resume-partial.
if ((resume_partial)); then
    submit_metals=("${metals[@]}")
else
    submit_metals=("${created[@]}")
    if ((include_existing)); then
        for metal in "${metals[@]}"; do
            for newmetal in "${created[@]}"; do
                [[ "$metal" != "$newmetal" ]] || continue 2
            done
            submit_metals+=("$metal")
        done
    fi
fi

submit_metal() (
    local metal="$1" dir="$base/$1" stage jid dep='' memory time_limit raw cpus partition
    local start_index=0 result status recorded_count saved_dep i
    if [[ ! -f "$dir/metal.json" || ! -f "$dir/control" || ! -f "$dir/coord" || \
          ! -f "$dir/def.inp" || ! -f "$dir/tm_workflow.py" || \
          ! -f "$dir/turbo.pipeline.optimized.slurm" || -e "$dir/control.template" ]]; then
        echo "ERROR SUBMIT $metal: not a valid prepared folder, or control.template exists" >&2
        return 1
    fi
    # Never submit consolidated stage names to an older per-metal script.
    if ! grep -Fq 'ridft_eiger_fix' "$dir/turbo.pipeline.optimized.slurm" ||
       ! grep -Fq 'jobex_eiger' "$dir/turbo.pipeline.optimized.slurm"; then
        echo "ERROR SUBMIT $metal: per-metal SLURM script is outdated; review/update it manually." >&2
        return 1
    fi

    # Serializes two simultaneous lan-run processes targeting the same metal.
    # A stale lock after SIGKILL requires inspection; never break it automatically.
    if ! mkdir "$dir/.lan_submission_lock" 2>/dev/null; then
        echo "ERROR SUBMIT $metal: another submission may be in progress (.lan_submission_lock exists)" >&2
        return 1
    fi
    trap 'rmdir "$dir/.lan_submission_lock" 2>/dev/null || true' EXIT
    if [[ -e "$dir/.lan_submission_pending" ]]; then
        echo "ERROR SUBMIT $metal: unresolved .lan_submission_pending; check scheduler before proceeding" >&2
        return 1
    fi

    local stages=(ridft_eiger_fix integer_ridft jobex_eiger)
    if [[ "$frequency" != none ]]; then stages+=("$frequency"); fi

    if ((resume_partial)); then
        if ! result="$(python3 "$SCRIPT_DIR/submission_guard.py" resume --dir "$dir" --stages "${stages[@]}")"; then
            echo "ERROR RESUME $metal: previous submission cannot safely be resumed" >&2
            return 1
        fi
        if [[ "$result" == FULL ]]; then
            echo "SKIP RESUME $metal: all ${#stages[@]} stages were already submitted"
            return 0
        fi
        IFS=$'\t' read -r status recorded_count saved_dep <<< "$result"
        if [[ "$status" != PARTIAL || ! "$recorded_count" =~ ^[0-9]+$ ]]; then
            echo "ERROR RESUME $metal: unrecognized verification response: $result" >&2
            return 1
        fi
        start_index="$recorded_count"
        if [[ "$saved_dep" != '-' ]]; then dep="$saved_dep"; fi
        echo "RESUME $metal: ${start_index}/${#stages[@]} stages already submitted; submitting only missing stages"
    else
        if [[ -e "$dir/chain_submission.tsv" || -e "$dir/chain_plan.txt" ]]; then
            echo "SKIP SUBMIT $metal: chain record exists (use --resume-partial only if incomplete)"
            return 0
        fi
        # This also protects against edited control, previously generated orbitals,
        # output files, and any unknown pre-existing contents; does not edit them.
        local prepared_opts=(prepared --base "$base" --metal "$metal")
        # Only already-existing folders need this extra check; newly created
        # ones cannot have pending pre-existing jobs by construction.
        local is_new=0 newmetal
        for newmetal in "${created[@]}"; do
            [[ "$metal" != "$newmetal" ]] || { is_new=1; break; }
        done
        if (( ! is_new )); then prepared_opts+=(--check-queue); fi
        if ! python3 "$SCRIPT_DIR/submission_guard.py" "${prepared_opts[@]}" >/dev/null; then
            echo "ERROR SUBMIT $metal: existing inputs are not demonstrably pristine; no files changed" >&2
            return 1
        fi
    fi

    # Each successful sbatch is journaled immediately. The pending marker
    # guards the unavoidable crash window between sbatch accepting a job and
    # writing its ID to disk. If the result is ambiguous, refuse blind replay.
    for ((i=start_index; i<${#stages[@]}; i++)); do
        stage="${stages[i]}"
        partition='general'
        case "$stage" in
            ridft_eiger_fix|integer_ridft) time_limit='2-00:00:00'; memory="$RIDFT_MEM"; cpus="$TM_CPUS" ;;
            jobex_eiger) time_limit='2-00:00:00'; memory="$JOBEX_MEM"; cpus="$TM_CPUS" ;;
            numforce|aoforce) time_limit='2-00:00:00'; memory="$HESS_MEM"; cpus="$TM_CPUS" ;;
        esac
        local export_spec="ALL,TM_STAGE=$stage"
        if [[ "$frequency" == none && "$stage" == jobex_eiger ]]; then
            export_spec+=",TM_FINAL_NO_HESSIAN=1"
        fi
        local opts=(--parsable --partition="$partition" --time="$time_limit"
                    --mem="$memory" --cpus-per-task="$cpus"
                    --job-name="tm-${metal}-${stage}"
                    --export="$export_spec")
        if [[ -n "$dep" ]]; then opts+=("--dependency=afterok:$dep"); fi
        printf 'stage=%s\nmetal=%s\n' "$stage" "$metal" > "$dir/.lan_submission_pending" || return 1
        if ! raw="$(cd "$dir" && sbatch "${opts[@]}" turbo.pipeline.optimized.slurm)"; then
            rm -f "$dir/.lan_submission_pending"
            echo "ERROR $metal: sbatch failed at $stage; previous submitted stages are preserved." >&2
            return 1
        fi
        jid="${raw%%;*}"
        if [[ ! "$jid" =~ ^[0-9]+$ ]]; then
            echo "ERROR $metal: unexpected sbatch output for $stage: $raw; pending marker retained for review" >&2
            return 1
        fi
        if [[ "$i" == 0 ]]; then
            # Write the selected full plan once. Prevents silent mode changes
            # when resuming an incomplete chain (especially optional Hessian).
            printf '%s\n' "${stages[@]}" > "$dir/chain_plan.txt.tmp.$$" || return 1
            mv "$dir/chain_plan.txt.tmp.$$" "$dir/chain_plan.txt" || return 1
            printf 'stage\tjob_id\tpartition\twalltime\tmemory\n' > "$dir/chain_submission.tsv" || return 1
        fi
        printf '%s\t%s\t%s\t%s\t%s\n' "$stage" "$jid" "$partition" "$time_limit" "$memory" >> "$dir/chain_submission.tsv" || return 1
        rm -f "$dir/.lan_submission_pending" || return 1
        printf '%-2s %-14s %-10s %s\n' "$metal" "$stage" "$jid" "${dep:+afterok:$dep}"
        dep="$jid"
    done
    if [[ "$frequency" == none ]]; then
        echo "SUBMITTED $metal: 3-job geometry-only chain (6 scientific stages; no Hessian; final SLURM job $dep)"
    else
        echo "SUBMITTED $metal: 4-job $frequency chain (7 scientific stages; final SLURM job $dep)"
    fi
)

for metal in "${submit_metals[@]}"; do
    if ! submit_metal "$metal"; then
        ((failures+=1))
        echo "WARNING: continuing after $metal submission failure" >&2
    fi
done

printf 'Done. Prepared %d new folder(s); %d failure(s). See each metal/chain_submission.tsv.\n' "${#created[@]}" "$failures"
echo 'To monitor: squeue -u "$USER"'
echo 'On failure, afterok blocks all dependent stages; orbital holes are warnings only.'
if ((failures > 0)); then exit 1; fi
