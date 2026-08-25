#!/usr/bin/env bash
#SBATCH --job-name=step17_archive
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=/path/to/mollusc_project/logs/step17_archive_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step17_archive_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
: "${STEP17_JOB_ID:?STEP17_JOB_ID is required}"
SOURCE=$PROJECT/results_or_reports/reproducibility_result_archiving/final_v1/run_${STEP17_JOB_ID}
test -f "$SOURCE/PROJECT_ARCHIVE_COMPLETE"
test -f "$SOURCE/COMPLETE"
"$PROJECT/scripts/28_sync_results_archive.sh" one \
  step17_reproducibility_result_archiving "$SOURCE" "run_${STEP17_JOB_ID}"

