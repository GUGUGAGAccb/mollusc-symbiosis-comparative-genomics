#!/usr/bin/env bash
#SBATCH --job-name=step11_combine
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step11_combine_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step11_combine_%j.err

set -euo pipefail

ROOT=/path/to/mollusc_project
BIV="$ROOT/results_or_reports/gene_tree_species_tree_reconciliation/bivalvia/dl_lca_v1/latest"
GAS="$ROOT/results_or_reports/gene_tree_species_tree_reconciliation/gastropoda/dl_lca_v1/latest"
BASE="$ROOT/results_or_reports/gene_tree_species_tree_reconciliation/combined_summary/dl_lca_v1"
RUN="$BASE/run_${SLURM_JOB_ID}"
SCRIPT="$ROOT/scripts/27_step11_reconciliation.py"
PYTHON="$ROOT/conda_envs/orthofinder_env/bin/python"

mkdir -p "$RUN" "$ROOT/logs"
test -s "$BIV/COMPLETE"
test -s "$GAS/COMPLETE"

"$PYTHON" "$SCRIPT" combine \
  --bivalvia "$BIV" \
  --gastropoda "$GAS" \
  --output-dir "$RUN"

ln -sfn "run_${SLURM_JOB_ID}" "$BASE/latest"
ARCHIVER="$ROOT/scripts/28_sync_results_archive.sh"
if [[ -f "$ARCHIVER" ]]; then
  bash "$ARCHIVER" one "step11_gene_tree_species_tree_reconciliation/bivalvia" "$BIV" "$(basename "$(readlink -f "$BIV")")"
  bash "$ARCHIVER" one "step11_gene_tree_species_tree_reconciliation/gastropoda" "$GAS" "$(basename "$(readlink -f "$GAS")")"
  bash "$ARCHIVER" one "step11_gene_tree_species_tree_reconciliation/combined" "$RUN" "$(basename "$RUN")"
fi
printf 'Step 11 combined summary complete: %s\n' "$RUN"
