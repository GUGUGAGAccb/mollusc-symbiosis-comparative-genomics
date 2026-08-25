#!/bin/bash
#SBATCH --job-name=orthogroup_distribution
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step04_orthogroup_distribution_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step04_orthogroup_distribution_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
SCRIPT="$PROJECT/scripts/16_step04_orthogroup_distribution.py"
OUTPUT_ROOT="$PROJECT/results_or_reports/orthogroup_distribution"
RUN_DIR="$OUTPUT_ROOT/run_${SLURM_JOB_ID}"

mkdir -p "$OUTPUT_ROOT" "$PROJECT/logs"
python3 "$SCRIPT" --project-root "$PROJECT" --output "$RUN_DIR"
ln -sfn "$RUN_DIR" "$OUTPUT_ROOT/latest"

test -s "$RUN_DIR/COMPLETE"
test -s "$RUN_DIR/step04_summary.tsv"
echo "STEP04_COMPLETE=$RUN_DIR"
