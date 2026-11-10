#!/bin/bash
# VACC/TURBOMOLE 7.9: one SLURM allocation for ONE calculation stage.
# Do not sbatch this file directly: submit-turbo-chain.optimized.sh supplies
# TM_STAGE and creates the afterok dependency chain.
#SBATCH --job-name=tm-stage
#SBATCH --account=pi-crreinha
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --mail-type=NONE

# Load Lmod before enabling nounset: VACC module initialization may read
# unset shell variables (notably PS1 in noninteractive SLURM jobs).
set -Eeo pipefail
module purge
module load turbomole/7.9
if [[ -z "${TURBODIR:-}" ]]; then
    echo "ERROR: turbomole/7.9 was loaded, but TURBODIR is unset; inspect 'module show turbomole/7.9'." >&2
    exit 1
fi
set -u
ulimit -s unlimited

# VERIFIED VACC NumForce/OpenMP settings -- do not change this configuration.
unset PARA_ARCH TM_PAR_OMP TM_PAR_FORK SMP SMPCPUS HOSTS_FILE
export PATH="$TURBODIR/bin/em64t-unknown-linux-gnu:$TURBODIR/scripts:$PATH"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export PARNODES="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores
export OMP_PROC_BIND=close

export TURBOTMPDIR="$HOME/scratch/turbomole_tmp/$SLURM_JOB_ID"
mkdir -p "$TURBOTMPDIR"
cd "$SLURM_SUBMIT_DIR"

[[ -f control && -f coord && -f tm_workflow.py && -f metal.json ]] || {
    echo "ERROR: control, coord, tm_workflow.py or metal.json missing in $PWD" >&2
    exit 1
}
[[ ! -e control.template ]] || {
    echo "ERROR: control.template exists; DEFINE needs the single control file" >&2
    exit 1
}

stage="${TM_STAGE:-}"
case "$stage" in
    ridft|eiger_pre|occupation_fix|integer_ridft|jobex|eiger_post|numforce|aoforce) ;;
    *) echo "ERROR: unknown TM_STAGE='$stage'; use the chain submission script." >&2; exit 2 ;;
esac

on_exit() {
    local rc=$?
    if (( rc != 0 )); then
        python3 tm_workflow.py record "${stage^^}" FAILED "SLURM stage exited with code $rc; see slurm output" || true
        echo "ERROR: ${stage} failed (exit $rc)" >&2
    fi
}
trap on_exit EXIT

printf '\n======= TURBOMOLE STAGE =======\n'
date
printf 'Stage: %s\nJob: %s\nNode: %s\nPartition: %s\nCPUs: %s\nOMP: %s\nPARNODES: %s\nDirectory: %s\nScratch: %s\n' \
    "$stage" "$SLURM_JOB_ID" "$SLURM_JOB_NODELIST" \
    "$SLURM_JOB_PARTITION" "$SLURM_CPUS_PER_TASK" \
    "$OMP_NUM_THREADS" "$PARNODES" "$PWD" "$TURBOTMPDIR"

case "$stage" in
    ridft)
        # The sole control file exists BEFORE define. No control.template,
        # no define.out and no rewriting control after define.
        [[ -f def.inp ]] || { echo 'ERROR: def.inp missing' >&2; exit 1; }
        python3 tm_workflow.py record DEFINE RUNNING 'define < def.inp'
        define < def.inp
        python3 tm_workflow.py check-define

        actual -r
        python3 tm_workflow.py record RIDFT RUNNING 'initial SCF with Fermi smearing'
        if ! ridft > ridft.out 2>&1 || ! python3 tm_workflow.py check-ridft ridft.out; then
            python3 tm_workflow.py record RIDFT RETRY 'one retry with the existing orbitals'
            if ! ridft > ridft.retry.out 2>&1 || ! python3 tm_workflow.py check-ridft ridft.retry.out; then
                echo 'ERROR: RIDFT did not converge in the initial run or retry.' >&2
                exit 1
            fi
        fi
        ;;

    eiger_pre)
        # An occupation warning or failed eiger does NOT stop geometry optimization.
        if ! eiger > eiger.out 2>&1; then
            echo 'WARNING: eiger failed before geometry optimization; continuing' >&2
        fi
        python3 tm_workflow.py eiger eiger.out pre || true
        ;;

    occupation_fix)
        # Quick input edit only: restore initial alpha/beta shell assignments,
        # retain the previous workflow's removal of $fermi, and write control.
        # No TURBOMOLE SCF is run in this short, 1-CPU, 1-GB allocation.
        python3 tm_workflow.py remove-fermi
        ;;

    integer_ridft)
        # Run the same integer-occupation SCF as previously, now separately
        # on general so that it retains up to 48 hours and 20 CPUs.
        ridft > ridft.integer.out 2>&1
        python3 tm_workflow.py check-ridft ridft.integer.out
        python3 tm_workflow.py restore-occupations
        ;;

    jobex)
        # Clear stale success markers only; preserve the actual molecular data.
        rm -f GEO_OPT_CONVERGED GEO_OPT_FAILED converged not.converged
        python3 tm_workflow.py record JOBEX RUNNING 'jobex -ri -gcart 4 -c 600'
        jobex -ri -gcart 4 -c 600 > jobex.out 2>&1
        if [[ -e GEO_OPT_FAILED || -e not.converged ]] || \
           { [[ ! -f GEO_OPT_CONVERGED ]] && [[ ! -f converged ]]; }; then
            echo 'ERROR: JOBEX did not reach geometry convergence.' >&2
            exit 2
        fi
        # Standard TURBOMOLE may write 'converged' rather than the earlier
        # package's GEO_OPT_CONVERGED marker. Normalize ONLY after success.
        [[ -f GEO_OPT_CONVERGED ]] || touch GEO_OPT_CONVERGED
        python3 tm_workflow.py record JOBEX OK 'geometry converged'
        ;;

    eiger_post)
        if ! eiger > eiger.after_jobex.out 2>&1; then
            echo 'WARNING: eiger failed after jobex; continuing' >&2
        fi
        python3 tm_workflow.py eiger eiger.after_jobex.out post || true
        ;;

    numforce)
        [[ -f GEO_OPT_CONVERGED && -s gradient ]] || {
            echo 'ERROR: successful optimization marker or gradient missing.' >&2
            exit 2
        }
        python3 tm_workflow.py record NUMFORCE RUNNING 'PARNODES=1 NumForce -ri -central'
        # VERIFIED: exactly one displacement worker, OpenMP threads from allocation.
        PARNODES=1 NumForce -ri -central > numforce.out 2>&1
        python3 tm_workflow.py check-numforce
        python3 tm_workflow.py record DONE SUCCESS 'NumForce completed; see occupation warnings'
        ;;

    aoforce)
        [[ -f GEO_OPT_CONVERGED ]] || {
            echo 'ERROR: successful optimization marker missing.' >&2
            exit 2
        }
        python3 tm_workflow.py record AOFORCE RUNNING 'aoforce'
        aoforce > aoforce.out 2>&1
        python3 tm_workflow.py check-aoforce
        python3 tm_workflow.py record DONE SUCCESS 'aoforce completed; see occupation warnings'
        ;;
esac

echo "Completed stage $stage at $(date)"
