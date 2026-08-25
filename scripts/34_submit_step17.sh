#!/usr/bin/env bash
set -euo pipefail

PROJECT=/path/to/mollusc_project
main_job=$(sbatch --parsable "$PROJECT/scripts/34_step17_reproducibility_archive.slurm")
archive_job=$(sbatch --parsable --dependency=afterok:"$main_job" --export=ALL,STEP17_JOB_ID="$main_job" "$PROJECT/scripts/34_step17_archive.slurm")
printf '%s\t%s\n' "$main_job" "$archive_job"

