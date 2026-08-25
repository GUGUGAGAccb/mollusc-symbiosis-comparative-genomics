#!/usr/bin/env bash
#SBATCH --job-name=step14_archive
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step14_archive_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step14_archive_%j.err

set -euo pipefail
PRIORITY_JOB_ID=${1:?Usage: sbatch script PRIORITY_JOB_ID}
ROOT=/path/to/mollusc_project
RUN="$ROOT/results_or_reports/integrated_candidate_prioritisation/evidence_weighted_v1/run_${PRIORITY_JOB_ID}"
test -s "$RUN/COMPLETE"
bash "$ROOT/scripts/28_sync_results_archive.sh" one \
  step14_integrated_candidate_prioritisation "$RUN" "run_${PRIORITY_JOB_ID}"
