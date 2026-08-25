#!/usr/bin/env bash
#SBATCH --job-name=step16_stats_viz
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/path/to/mollusc_project/logs/step16_stats_viz_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step16_stats_viz_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
STEP14_RUN=${STEP14_RUN:-$PROJECT/results_or_reports/integrated_candidate_prioritisation/evidence_weighted_v1/run_18427340}
STEP15_RUN=${STEP15_RUN:-$PROJECT/results_or_reports/sensitivity_bias_analysis/ranking_sensitivity_v1/run_18431248}
OUT_ROOT=$PROJECT/results_or_reports/statistical_analysis_visualisation/integrated_figures_v1
OUT=$OUT_ROOT/run_${SLURM_JOB_ID}
RSCRIPT=$PROJECT/conda_envs/busco_env/bin/Rscript

test -f "$STEP14_RUN/COMPLETE"
test -f "$STEP15_RUN/COMPLETE"
test -x "$RSCRIPT"
mkdir -p "$OUT_ROOT"
"$RSCRIPT" "$PROJECT/scripts/33_step16_statistics_visualisation.R" \
  "$PROJECT" "$STEP14_RUN" "$STEP15_RUN" "$OUT"
ln -sfn "run_${SLURM_JOB_ID}" "$OUT_ROOT/latest"

