#!/usr/bin/env bash
set -euo pipefail
ROOT=/path/to/mollusc_project
priority=$(sbatch --parsable "$ROOT/scripts/31_step14_prioritise.slurm")
archive=$(sbatch --parsable --dependency=afterok:"$priority" \
  "$ROOT/scripts/31_step14_archive.slurm" "$priority")
printf 'role\tjob_id\nprioritisation\t%s\narchive\t%s\n' "$priority" "$archive"
