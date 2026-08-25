#!/usr/bin/env bash
set -euo pipefail

PROJECT=/path/to/mollusc_project
STEP14_RUN=${1:-$PROJECT/results_or_reports/integrated_candidate_prioritisation/evidence_weighted_v1/run_18427340}

test -f "$STEP14_RUN/COMPLETE"
main_job=$(sbatch --parsable --export=ALL,STEP14_RUN="$STEP14_RUN" "$PROJECT/scripts/32_step15_sensitivity.slurm")
archive_job=$(sbatch --parsable --dependency=afterok:"$main_job" --export=ALL,STEP15_JOB_ID="$main_job" "$PROJECT/scripts/32_step15_archive.slurm")
printf '%s\t%s\n' "$main_job" "$archive_job"
