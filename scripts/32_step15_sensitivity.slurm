#!/usr/bin/env bash
#SBATCH --job-name=step15_sensitivity
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/path/to/mollusc_project/logs/step15_sensitivity_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step15_sensitivity_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
STEP14_RUN=${STEP14_RUN:-$PROJECT/results_or_reports/integrated_candidate_prioritisation/evidence_weighted_v1/run_18427340}
OUT_ROOT=$PROJECT/results_or_reports/sensitivity_bias_analysis/ranking_sensitivity_v1
OUT=$OUT_ROOT/run_${SLURM_JOB_ID}

mkdir -p "$OUT_ROOT"
python3 "$PROJECT/scripts/32_step15_sensitivity_bias.py" \
  --step14-run "$STEP14_RUN" \
  --busco "$PROJECT/primary_proteome_busco_summary.tsv" \
  --out "$OUT"
ln -sfn "run_${SLURM_JOB_ID}" "$OUT_ROOT/latest"
