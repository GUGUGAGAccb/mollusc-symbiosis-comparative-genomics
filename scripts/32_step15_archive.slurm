#!/usr/bin/env bash
#SBATCH --job-name=step15_archive
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/path/to/mollusc_project/logs/step15_archive_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step15_archive_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
: "${STEP15_JOB_ID:?STEP15_JOB_ID is required}"
SOURCE=$PROJECT/results_or_reports/sensitivity_bias_analysis/ranking_sensitivity_v1/run_${STEP15_JOB_ID}

test -f "$SOURCE/COMPLETE"
"$PROJECT/scripts/28_sync_results_archive.sh" one \
  step15_sensitivity_bias_analysis "$SOURCE" "run_${STEP15_JOB_ID}"
