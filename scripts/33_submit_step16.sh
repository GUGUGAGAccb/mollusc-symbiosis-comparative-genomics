#!/usr/bin/env bash
set -euo pipefail

PROJECT=/path/to/mollusc_project
STEP14_RUN=${1:-$PROJECT/results_or_reports/integrated_candidate_prioritisation/evidence_weighted_v1/run_18427340}
STEP15_RUN=${2:-$PROJECT/results_or_reports/sensitivity_bias_analysis/ranking_sensitivity_v1/run_18431248}
test -f "$STEP14_RUN/COMPLETE"
test -f "$STEP15_RUN/COMPLETE"
main_job=$(sbatch --parsable --export=ALL,STEP14_RUN="$STEP14_RUN",STEP15_RUN="$STEP15_RUN" "$PROJECT/scripts/33_step16_statistics_visualisation.slurm")
archive_job=$(sbatch --parsable --dependency=afterok:"$main_job" --export=ALL,STEP16_JOB_ID="$main_job" "$PROJECT/scripts/33_step16_archive.slurm")
printf '%s\t%s\n' "$main_job" "$archive_job"

