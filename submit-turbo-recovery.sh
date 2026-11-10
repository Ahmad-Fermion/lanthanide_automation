#!/bin/bash
# Recovery submitter: restart RIDFT after manual edits, optionally link remaining stages.
# Does not call define and never creates or modifies input folders.
set -Eeuo pipefail

usage() {
cat <<'HELP'
Restart RIDFT after manual SCF adjustments in an existing metal folder.
No define or input replacement. Run from the EF project directory.

Usage: lan-recover [num|ao] --metal METAL [mode]   (default: num)
  num = NumForce; ao = aoforce (full names also accepted).

Modes (choose one):
  --ridft-only     Restart only RIDFT (default); inspect results manually
  --continue-only  After a successful recovery, submit remaining stages
  --full-chain     Restart RIDFT and schedule all remaining dependent stages

Other options: --base DIR (default: current directory); -h, --help

Examples:
  lan-recover num --metal tb                  # RIDFT only
  lan-recover num --metal tb --continue-only  # continue after successful RIDFT
  lan-recover ao --metal gd --full-chain      # RIDFT + all following stages

Check/cancel old queued jobs for this metal before recovery. Other metals
are unaffected. Later occupation_fix still restores occupations/removes $fermi.
HELP
}

mode='ridft-only'
frequency='numforce'
base="$PWD"
metal=''
while (($#)); do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        num|numforce) frequency='numforce'; shift ;;
        ao|aoforce) frequency='aoforce'; shift ;;
        --metal)
            (($# >= 2)) || { echo 'ERROR: --metal needs a symbol' >&2; exit 2; }
            metal="${2,,}"; shift 2 ;;
        --base)
            (($# >= 2)) || { echo 'ERROR: --base needs a directory' >&2; exit 2; }
            base="$2"; shift 2 ;;
        --ridft-only) mode='ridft-only'; shift ;;
        --full-chain) mode='full-chain'; shift ;;
        --continue-only) mode='continue-only'; shift ;;
        *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done
case "$metal" in
    la|ce|pr|nd|sm|eu|gd|tb|dy|ho|er|tm|yb|lu) ;;
    *) echo 'ERROR: specify exactly one supported metal with --metal (Pm excluded)' >&2; exit 2 ;;
esac
base="$(cd -- "$base" && pwd -P)"
dir="$base/$metal"
[[ -d "$dir" ]] || { echo "ERROR: $dir is missing" >&2; exit 1; }
for name in control coord metal.json tm_workflow.py turbo.pipeline.optimized.slurm; do
    [[ -f "$dir/$name" ]] || { echo "ERROR: missing $dir/$name" >&2; exit 1; }
done
if [[ "$mode" != continue-only ]]; then
    [[ -f "$dir/turbo.resume.slurm" ]] || {
        echo "ERROR: missing $dir/turbo.resume.slurm" >&2
        echo "Copy it with: cp lanthanide_automation/turbo.resume.slurm $metal/" >&2
        exit 1
    }
fi
if [[ "$mode" == continue-only ]]; then
    [[ -f "$dir/ridft.resume.ok" ]] || {
        echo 'ERROR: ridft.resume.ok missing. Finish and verify RIDFT recovery first.' >&2; exit 1;
    }
fi
if [[ "$mode" != ridft-only ]]; then
    [[ -s "$dir/initial_occupations" ]] || {
        echo 'ERROR: initial_occupations missing; the normal occupation_fix stage requires it.' >&2; exit 1;
    }
fi

# Check the SLURM WORK DIRECTORY, not the job name. The user may have
# identically named Sm jobs in several different complexes/models at once.
if command -v squeue >/dev/null 2>&1; then
    jobs="$(squeue -h -u "${USER:-}" -o '%i|%Z')" || {
        echo 'ERROR: could not verify SLURM queue state; refusing ambiguous submission.' >&2
        exit 1
    }
    while IFS='|' read -r queued_id queued_workdir; do
        [[ -n "$queued_id" ]] || continue
        if [[ "$queued_workdir" == "$dir" ]]; then
            echo "ERROR: job $queued_id is still queued/running in $dir." >&2
            echo 'Inspect/cancel obsolete dependency jobs before submitting a recovery.' >&2
            exit 1
        fi
    done <<< "$jobs"
fi

: "${TM_CPUS:=20}" "${EIGER_CPUS:=1}" "${RIDFT_MEM:=20G}" "${JOBEX_MEM:=20G}" "${HESS_MEM:=100G}" "${EIGER_MEM:=4G}"
: "${RECOVERY_CPUS:=$TM_CPUS}" "${RECOVERY_MEM:=$RIDFT_MEM}" "${RECOVERY_TIME:=2-00:00:00}"
for var in TM_CPUS EIGER_CPUS RECOVERY_CPUS; do
    [[ "${!var}" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid $var=${!var}" >&2; exit 2; }
done

record="$dir/recovery_submission.$(date -u +%Y%m%dT%H%M%SZ).$$.tsv"
printf 'stage\tjob_id\tpartition\twalltime\tmemory\tdependency\n' > "$record"
dep=''
submit_one() {
    local stage="$1" script="$2" partition="$3" time="$4" mem="$5" cpus="$6"
    local raw jid
    local args=(--parsable --job-name="tm-${metal}-${stage}" --partition="$partition"
                --time="$time" --mem="$mem" --cpus-per-task="$cpus")
    if [[ "$script" == 'turbo.pipeline.optimized.slurm' ]]; then
        args+=(--export="ALL,TM_STAGE=$stage")
    fi
    if [[ -n "$dep" ]]; then args+=(--dependency="afterok:$dep"); fi
    raw="$(cd "$dir" && sbatch "${args[@]}" "$script")" || {
        echo "ERROR: failed to submit $stage; see partial record $record" >&2; return 1;
    }
    jid="${raw%%;*}"
    [[ "$jid" =~ ^[0-9]+$ ]] || { echo "ERROR: invalid sbatch output: $raw" >&2; return 1; }
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$stage" "$jid" "$partition" "$time" "$mem" "${dep:+afterok:$dep}" >> "$record"
    printf '%-16s %s %s\n' "$stage" "$jid" "${dep:+afterok:$dep}"
    dep="$jid"
}

if [[ "$mode" != 'continue-only' ]]; then
    submit_one ridft-recovery turbo.resume.slurm general "$RECOVERY_TIME" "$RECOVERY_MEM" "$RECOVERY_CPUS"
fi
if [[ "$mode" != 'ridft-only' ]]; then
    for stage in eiger_pre occupation_fix integer_ridft jobex eiger_post "$frequency"; do
        partition=general
        case "$stage" in
            eiger_*) time='00:10:00'; mem="$EIGER_MEM"; cpus="$EIGER_CPUS" ;;
            occupation_fix) partition=short; time='00:10:00'; mem=1G; cpus=1 ;;
            integer_ridft) time='2-00:00:00'; mem="$RIDFT_MEM"; cpus="$TM_CPUS" ;;
            jobex) time='2-00:00:00'; mem="$JOBEX_MEM"; cpus="$TM_CPUS" ;;
            numforce|aoforce) time='2-00:00:00'; mem="$HESS_MEM"; cpus="$TM_CPUS" ;;
        esac
        submit_one "$stage" turbo.pipeline.optimized.slurm "$partition" "$time" "$mem" "$cpus"
    done
fi
printf 'Submitted %s recovery for %s; see %s\n' "$mode" "$metal" "$record"
echo 'No define command was submitted. Manual control edits are preserved through RIDFT recovery.'
