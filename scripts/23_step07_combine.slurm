#!/usr/bin/env bash
#SBATCH --job-name=step07_combine
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step07_combine_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step07_combine_%j.err

set -euo pipefail

ROOT=/path/to/mollusc_project
SCRIPT="$ROOT/scripts/23_step07_functional_annotation.py"
BASE="$ROOT/results_or_reports/functional_annotation"
RUN="$BASE/combined_summary/run_${SLURM_JOB_ID}"
mkdir -p "$RUN"

python "$SCRIPT" combine \
  --bivalvia "$BASE/bivalvia/eggnog_mapper_v2_1_12/latest" \
  --gastropoda "$BASE/gastropoda/eggnog_mapper_v2_1_12/latest" \
  --output-dir "$RUN"

ln -sfn "run_${SLURM_JOB_ID}" "$BASE/combined_summary/latest"
printf 'Step 7 combined summary complete: %s\n' "$RUN"
