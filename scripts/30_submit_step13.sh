#!/usr/bin/env bash
set -euo pipefail
ROOT=/path/to/mollusc_project
SCRIPTS="$ROOT/scripts"

prep=$(sbatch --parsable "$SCRIPTS/30_step13_prepare.slurm")
array=$(sbatch --parsable --dependency=afterok:"$prep" --array=0-5%2 \
  "$SCRIPTS/30_step13_salmon_array.slurm" "$prep")
summary=$(sbatch --parsable --dependency=afterok:"$array" \
  "$SCRIPTS/30_step13_summarise.slurm" "$prep")
archive=$(sbatch --parsable --dependency=afterok:"$summary" \
  "$SCRIPTS/30_step13_archive.slurm" "$prep")

printf 'role\tjob_id\nprepare\t%s\nsalmon_array\t%s\nsummary\t%s\narchive\t%s\n' \
  "$prep" "$array" "$summary" "$archive"
