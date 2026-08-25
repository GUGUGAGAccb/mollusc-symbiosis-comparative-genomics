#!/usr/bin/env bash
#SBATCH --job-name=step09_combine
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step09_combine_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step09_combine_%j.err

set -euo pipefail
ROOT=/path/to/mollusc_project
SCRIPT="$ROOT/scripts/25_step09_candidate_sequence_validation.py"
BASE="$ROOT/results_or_reports/candidate_sequence_validation"
RUN="$BASE/combined_summary/run_${SLURM_JOB_ID}"
mkdir -p "$RUN"
python "$SCRIPT" combine \
  --bivalvia "$BASE/bivalvia/pfam_hmmer_v1/latest" \
  --gastropoda "$BASE/gastropoda/pfam_hmmer_v1/latest" \
  --output-dir "$RUN"
ln -sfn "run_${SLURM_JOB_ID}" "$BASE/combined_summary/latest"
printf 'Step 9 combined summary complete: %s\n' "$RUN"
