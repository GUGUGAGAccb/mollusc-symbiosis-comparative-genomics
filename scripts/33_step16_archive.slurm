#!/usr/bin/env bash
#SBATCH --job-name=step16_archive
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/path/to/mollusc_project/logs/step16_archive_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step16_archive_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
: "${STEP16_JOB_ID:?STEP16_JOB_ID is required}"
SOURCE=$PROJECT/results_or_reports/statistical_analysis_visualisation/integrated_figures_v1/run_${STEP16_JOB_ID}
test -f "$SOURCE/COMPLETE"
"$PROJECT/scripts/28_sync_results_archive.sh" one \
  step16_statistical_analysis_visualisation "$SOURCE" "run_${STEP16_JOB_ID}"

