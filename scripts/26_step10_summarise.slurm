#!/usr/bin/env bash
#SBATCH --job-name=step10_summarise
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step10_summarise_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step10_summarise_%j.err

set -euo pipefail

ROOT=/path/to/mollusc_project
BASE="$ROOT/results_or_reports/candidate_family_phylogenetics/validated_top30_v1"
RUN=$(readlink -f "$BASE/latest_prepared")
SCRIPT="$ROOT/scripts/26_step10_candidate_gene_phylogenetics.py"

python "$SCRIPT" summarise --run-dir "$RUN"
ln -sfn "$(basename "$RUN")" "$BASE/latest"
printf 'Step 10 summary complete: %s\n' "$RUN"
