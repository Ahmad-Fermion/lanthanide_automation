#!/bin/bash
# VACC/TURBOMOLE 7.9: group lightweight stages into compute allocations.
# Do not sbatch this file directly: submit-turbo-chain.sh supplies
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

[[ -f control && -f coord && -f tm_workflow.py ]] &&
    { [[ -f metal.json ]] || [[ "${TM_MANUAL_RECOVERY:-0}" == 1 ]]; } || {
    echo "ERROR: control/coord/tm_workflow.py missing, or metal.json unavailable outside manual Hessian recovery: $PWD" >&2
    exit 1
}
[[ ! -e control.template ]] || {
    echo "ERROR: control.template exists; DEFINE needs the single control file" >&2
    exit 1
}

stage="${TM_STAGE:-}"
case "$stage" in
    ridft|ridft_eiger_fix|eiger_pre|occupation_fix|integer_ridft|eiger_fix_integer_ridft|jobex|jobex_eiger|screwer_jobex_eiger|screwer_ridft_jobex_eiger|restricted_ridft|restricted_ridft_no_define|eiger_post|numforce|aoforce) ;;
    *) echo "ERROR: unknown TM_STAGE='$stage'; use the chain submission script." >&2; exit 2 ;;
esac
if [[ "${TM_MANUAL_RECOVERY:-0}" == 1 && "$stage" != numforce && "$stage" != aoforce ]]; then
    echo 'ERROR: metadata-free manual recovery is limited to Hessian-only stages' >&2
    exit 2
fi

active_stage="$stage"
on_exit() {
    local rc=$?
    if (( rc != 0 )); then
        python3 tm_workflow.py record "${active_stage^^}" FAILED "SLURM stage exited with code $rc; see slurm output" || true
        echo "ERROR: ${active_stage} failed (exit $rc)" >&2
    fi
}
trap on_exit EXIT

printf '\n======= TURBOMOLE STAGE =======\n'
date
printf 'Stage: %s\nJob: %s\nNode: %s\nPartition: %s\nCPUs: %s\nOMP: %s\nPARNODES: %s\nDirectory: %s\nScratch: %s\n' \
    "$stage" "$SLURM_JOB_ID" "$SLURM_JOB_NODELIST" \
    "$SLURM_JOB_PARTITION" "$SLURM_CPUS_PER_TASK" \
    "$OMP_NUM_THREADS" "$PARNODES" "$PWD" "$TURBOTMPDIR"

# Preserve the original scientific stages while running quick checks/edits
# inside the adjacent compute allocation. EIGER warnings stay non-blocking.
run_eiger_pre() {
    active_stage=eiger_pre
    if ! eiger > eiger.out 2>&1; then
        echo 'WARNING: eiger failed before geometry optimization; continuing' >&2
    fi
    python3 tm_workflow.py eiger eiger.out pre || true
}

run_occupation_fix() {
    active_stage=occupation_fix
    python3 tm_workflow.py remove-fermi
}

run_eiger_post() {
    active_stage=eiger_post
    if ! eiger > eiger.after_jobex.out 2>&1; then
        echo 'WARNING: eiger failed after jobex; continuing' >&2
    fi
    python3 tm_workflow.py eiger eiger.after_jobex.out post || true
    if [[ "${TM_FINAL_NO_HESSIAN:-0}" == 1 ]]; then
        python3 tm_workflow.py record DONE SUCCESS 'Geometry optimization completed; Hessian intentionally skipped'
    fi
}

case "$stage" in
    ridft|ridft_eiger_fix)
        active_stage=ridft
        # The sole control file exists BEFORE define. No control.template,
        # no define.out and no rewriting control after define.
        define_input='def.inp'
        if [[ "${TM_SCREWER_CHILD:-0}" == 1 ]]; then
            [[ -f screwer_prepared.json ]] || { echo 'ERROR: screwer child manifest missing' >&2; exit 1; }
            define_input='def_screwer.inp'
        fi
        [[ -f "$define_input" ]] || { echo "ERROR: $define_input missing" >&2; exit 1; }
        python3 tm_workflow.py record DEFINE RUNNING "define < $define_input"
        define < "$define_input"
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
        if [[ "$stage" == ridft_eiger_fix ]]; then
            run_eiger_pre
            run_occupation_fix
        fi
        ;;

    restricted_ridft|restricted_ridft_no_define)
        # A distinct, opt-in RKS calculation: no Fermi or open-shell reset.
        [[ -s restricted_prepared.json ]] || { echo 'ERROR: restricted manifest missing' >&2; exit 2; }
        if [[ "$stage" == restricted_ridft ]]; then
            [[ -s def.inp ]] || { echo 'ERROR: restricted DEFINE input missing' >&2; exit 2; }
            active_stage=define
            python3 tm_workflow.py record DEFINE RUNNING 'restricted DEFINE'
            define < def.inp > define.out 2>&1
        fi
        python3 tm_workflow.py check-restricted-define
        active_stage=ridft
        actual -r
        python3 tm_workflow.py record RIDFT RUNNING 'restricted ground-state DFT'
        log=ridft.out
        [[ "$stage" == restricted_ridft ]] || log="ridft.resume.${SLURM_JOB_ID}.out"
        ridft > "$log" 2>&1
        python3 tm_workflow.py check-ridft "$log"
        if [[ "${TM_FINAL_NO_HESSIAN:-0}" == 1 ]]; then
            python3 tm_workflow.py record RIDFT OK 'restricted SCF finished; following JOBEX job will complete'
        fi
        ;;

    eiger_pre)
        run_eiger_pre
        ;;

    occupation_fix)
        run_occupation_fix
        ;;

    integer_ridft|eiger_fix_integer_ridft)
        if [[ "$stage" == eiger_fix_integer_ridft ]]; then
            run_eiger_pre
            run_occupation_fix
        fi
        active_stage=integer_ridft
        # Run the same integer-occupation SCF as previously, now separately
        # on general so that it retains up to 48 hours and 20 CPUs.
        ridft > ridft.integer.out 2>&1
        python3 tm_workflow.py check-ridft ridft.integer.out
        python3 tm_workflow.py restore-occupations
        ;;

    screwer_jobex_eiger|screwer_ridft_jobex_eiger)
        # Screwer continues with the converged parent orbitals: no Fermi
        # smearing, pre-JOBEX EIGER or occupation reset. Run one fresh SCF
        # for the displaced geometry before optimizing it.
        [[ "${TM_SCREWER_CHILD:-0}" == 1 && -f screwer_prepared.json && -f def_screwer.inp ]] || {
            echo 'ERROR: screwer continuation requires its manifest and def_screwer.inp' >&2
            exit 2
        }
        if [[ "$stage" == screwer_jobex_eiger ]]; then
        active_stage=define_screwer
        python3 tm_workflow.py record DEFINE RUNNING 'Screwer DEFINE with no-ired fallback on failure'
        python3 screwer_workflow.py --define-fallback
        python3 tm_workflow.py record DEFINE OK 'Screwer DEFINE validated; continuing without orbital/control comparison'
        else
            # Recovery must not execute Screwer displacement or DEFINE again.
            [[ -s control && -s coord ]] || { echo 'ERROR: missing inherited state' >&2; exit 2; }
        fi
        active_stage=ridft
        actual -r
        python3 tm_workflow.py record RIDFT RUNNING 'screwer continuation SCF before JOBEX'
        if [[ "$stage" == screwer_ridft_jobex_eiger ]]; then
            ridft > "ridft.resume.${SLURM_JOB_ID}.out" 2>&1
            python3 tm_workflow.py check-ridft "ridft.resume.${SLURM_JOB_ID}.out"
        else
            ridft > ridft.out 2>&1
            python3 tm_workflow.py check-ridft ridft.out
        fi
        # Fall through to the same JOBEX + EIGER implementation; RIDFT
        # shares the first SLURM allocation rather than requesting another.
        ;&

    jobex|jobex_eiger)
        active_stage=jobex
        # In recovery archive previous failure evidence before resetting markers.
        # In particular, not.converged MUST be deleted before restarted JOBEX.
        if [[ "${TM_RECOVERY:-0}" == 1 ]]; then
            archive="recovery_attempts/${SLURM_JOB_ID}"
            mkdir -p "$archive"
            for f in jobex.out GEO_OPT_FAILED not.converged GEO_OPT_CONVERGED converged; do
                [[ ! -f "$f" ]] || cp -p -- "$f" "$archive/$f"
            done
        fi
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
        if [[ "$stage" == jobex_eiger || "$stage" == screwer_jobex_eiger || "$stage" == screwer_ridft_jobex_eiger ]]; then
            run_eiger_post
        fi
        ;;

    eiger_post)
        run_eiger_post
        ;;

    numforce)
        [[ ( -f GEO_OPT_CONVERGED || -f converged ) && -s gradient ]] || {
            echo 'ERROR: successful optimization marker or gradient missing.' >&2
            exit 2
        }
        [[ -f GEO_OPT_CONVERGED ]] || touch GEO_OPT_CONVERGED
        if [[ "${TM_RECOVERY:-0}" == 1 ]]; then
            archive="recovery_attempts/${SLURM_JOB_ID}"
            mkdir -p "$archive"
            for f in hessian vibspectrum vib_normal_modes nprhessian nprvibspectrum npr_vib_normal_modes; do
                if [[ -f "$f" ]]; then mv -- "$f" "$archive/$f"; fi
            done
        elif [[ "${TM_SCREWER_CHILD:-0}" == 1 ]]; then
            python3 screwer_workflow.py --clean-frequency --base "$(dirname "$(dirname "$PWD")")" --metals "$(basename "$(dirname "$PWD")")"
        fi
        python3 tm_workflow.py record NUMFORCE RUNNING 'PARNODES=1 NumForce -ri -central'
        # VERIFIED: exactly one displacement worker, OpenMP threads from allocation.
        PARNODES=1 NumForce -ri -central > numforce.out 2>&1
        python3 tm_workflow.py check-numforce
        python3 tm_workflow.py record DONE SUCCESS 'NumForce completed; see occupation warnings'
        if [[ "${TM_SCREWER_CHILD:-0}" == 1 ]]; then
            python3 screwer_workflow.py --export-only --base "$(dirname "$(dirname "$PWD")")" --metals "$(basename "$(dirname "$PWD")")"
        fi
        ;;

    aoforce)
        [[ -f GEO_OPT_CONVERGED || -f converged ]] || {
            echo 'ERROR: successful optimization marker missing.' >&2
            exit 2
        }
        [[ -f GEO_OPT_CONVERGED ]] || touch GEO_OPT_CONVERGED
        if [[ "${TM_RECOVERY:-0}" == 1 ]]; then
            archive="recovery_attempts/${SLURM_JOB_ID}"
            mkdir -p "$archive"
            for f in hessian vibspectrum vib_normal_modes nprhessian nprvibspectrum npr_vib_normal_modes; do
                if [[ -f "$f" ]]; then mv -- "$f" "$archive/$f"; fi
            done
        elif [[ "${TM_SCREWER_CHILD:-0}" == 1 ]]; then
            python3 screwer_workflow.py --clean-frequency --base "$(dirname "$(dirname "$PWD")")" --metals "$(basename "$(dirname "$PWD")")"
        fi
        python3 tm_workflow.py record AOFORCE RUNNING 'aoforce'
        aoforce > aoforce.out 2>&1
        python3 tm_workflow.py check-aoforce
        python3 tm_workflow.py record DONE SUCCESS 'aoforce completed; see occupation warnings'
        if [[ "${TM_SCREWER_CHILD:-0}" == 1 ]]; then
            python3 screwer_workflow.py --export-only --base "$(dirname "$(dirname "$PWD")")" --metals "$(basename "$(dirname "$PWD")")"
        fi
        ;;
esac

echo "Completed stage $stage at $(date)"
