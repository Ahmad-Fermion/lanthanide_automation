#!/bin/bash
# Optional post-Hessian gate. Prepare one screwer distortion and submit
# its dependent optimization chain *only* if a verified imaginary mode exists.
#SBATCH --job-name=tm-screwer
#SBATCH --account=pi-crreinha
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --mail-type=NONE
set -Eeo pipefail
module purge
module load turbomole/7.9
if [[ -z "${TURBODIR:-}" ]]; then echo 'ERROR: TURBODIR not set' >&2; exit 1; fi
set -u
export PATH="$TURBODIR/bin/em64t-unknown-linux-gnu:$TURBODIR/scripts:$PATH"
cd "$SLURM_SUBMIT_DIR"
method="${TM_SCREWER_METHOD:?Missing Hessian method}"
metal="$(basename "$PWD")"
export OMP_NUM_THREADS=1
python3 screwer_workflow.py "$method" --base .. --metals "$metal" --min-imag "${TM_SCREWER_MIN_IMAG:-0}"
