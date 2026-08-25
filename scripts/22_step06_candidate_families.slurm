#!/bin/bash
#SBATCH --job-name=step06_candidates
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step06_candidates_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step06_candidates_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
SCRIPT="$PROJECT/scripts/22_step06_candidate_families.py"
OUTPUT_ROOT="$PROJECT/results_or_reports/candidate_gene_families/base_model_v1"
RUN_DIR="$OUTPUT_ROOT/run_${SLURM_JOB_ID}"

mkdir -p "$OUTPUT_ROOT" "$PROJECT/logs"
python3 "$SCRIPT" \
    --project "$PROJECT" \
    --output "$RUN_DIR" \
    --family-fdr 0.05 \
    --high-branch-p 0.01 \
    --candidate-branch-p 0.05 \
    --high-min-change 2

test -s "$RUN_DIR/COMPLETE"
test -s "$RUN_DIR/step06_summary.tsv"
ln -sfn "$RUN_DIR" "$OUTPUT_ROOT/latest"
echo "STEP06_COMPLETE=$RUN_DIR"
