#!/bin/bash
# MANUAL RIDFT RECOVERY ONLY: reuse the existing calculation and edited control.
# Do not run define, change $fermi, or restore reference occupations here.
# Submit from INSIDE a metal directory, e.g. sbatch turbo.resume.slurm.
#SBATCH --job-name=tm-ridft-recovery
#SBATCH --account=pi-crreinha
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --mail-type=NONE

# VACC Lmod can access unset PS1: only enable 'set -u' AFTER loading the module.
set -Eeo pipefail
module purge
module load turbomole/7.9
if [[ -z "${TURBODIR:-}" ]]; then
    echo 'ERROR: TURBODIR is unset after module load turbomole/7.9' >&2
    exit 1
fi
set -u
ulimit -s unlimited

# Do not alter tested VACC OpenMP setup.
unset PARA_ARCH TM_PAR_OMP TM_PAR_FORK SMP SMPCPUS HOSTS_FILE
export PATH="$TURBODIR/bin/em64t-unknown-linux-gnu:$TURBODIR/scripts:$PATH"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export PARNODES="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export TURBOTMPDIR="$HOME/scratch/turbomole_tmp/$SLURM_JOB_ID"
mkdir -p "$TURBOTMPDIR"
cd "$SLURM_SUBMIT_DIR"

for file in control coord metal.json tm_workflow.py; do
    [[ -f "$file" ]] || { echo "ERROR: $file missing in $PWD" >&2; exit 1; }
done
[[ ! -e control.template ]] || { echo 'ERROR: control.template must not exist' >&2; exit 1; }

# Every new run gets a distinct output. Keep prior SCF output and user edits.
log="ridft.resume.${SLURM_JOB_ID}.out"
rm -f ridft.resume.ok

on_exit() {
    local rc=$?
    if ((rc != 0)); then
        python3 tm_workflow.py record RIDFT_RECOVERY FAILED \
            "recovery job $SLURM_JOB_ID failed (exit $rc); see $log and SLURM log" || true
        echo "ERROR: RIDFT recovery failed; see $log" >&2
    fi
}
trap on_exit EXIT

python3 tm_workflow.py record RIDFT_RECOVERY RUNNING \
    "manual restart from edited control; output $log"
printf 'RIDFT recovery job=%s directory=%s threads=%s output=%s\n' \
    "$SLURM_JOB_ID" "$PWD" "$OMP_NUM_THREADS" "$log"
actual -r
ridft > "$log" 2>&1
python3 tm_workflow.py check-ridft "$log"
printf 'job_id=%s\nlog=%s\n' "$SLURM_JOB_ID" "$log" > ridft.resume.ok
python3 tm_workflow.py record RIDFT_RECOVERY OK \
    "SCF convergence verified; $log"
echo "RIDFT recovery converged: $log"
