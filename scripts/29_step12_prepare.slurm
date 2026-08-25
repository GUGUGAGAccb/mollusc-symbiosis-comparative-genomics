#!/usr/bin/env bash
#SBATCH --job-name=step12_prepare
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step12_prepare_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step12_prepare_%j.err

set -euo pipefail

ROOT=/path/to/mollusc_project
BASE="$ROOT/results_or_reports/molecular_evolution/hyphy_busted_absrel_relax_v1"
RUN="$BASE/run_${SLURM_JOB_ID}"
STEP10="$ROOT/results_or_reports/candidate_family_phylogenetics/validated_top30_v1/latest"
STEP11="$ROOT/results_or_reports/gene_tree_species_tree_reconciliation/combined_summary/dl_lca_v1/latest"
PYTHON="$ROOT/conda_envs/orthofinder_env/bin/python"
SCRIPT="$ROOT/scripts/29_step12_molecular_evolution.py"

mkdir -p "$RUN" "$ROOT/logs"
test -s "$STEP10/COMPLETE"
test -s "$STEP11/COMPLETE"
test -x "$PYTHON"
test -s "$SCRIPT"

"$PYTHON" "$SCRIPT" prepare \
  --root "$ROOT" \
  --step10 "$STEP10" \
  --step11 "$STEP11" \
  --output-dir "$RUN" \
  --min-sequences 5 \
  --min-species 3 \
  --min-codon-sites 100 \
  --min-mapping-fraction 0.70

test -s "$RUN/PREPARED"
ln -sfn "run_${SLURM_JOB_ID}" "$BASE/latest_prepared"
printf 'Step 12 preparation complete: %s\n' "$RUN"

