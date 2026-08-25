#!/usr/bin/env bash
set -euo pipefail

ROOT=/path/to/mollusc_project
SCRIPTS="$ROOT/scripts"

SETUP_JOB=$(sbatch --parsable "$SCRIPTS/29_step12_setup_hyphy.slurm")
PREP_JOB=$(sbatch --parsable "$SCRIPTS/29_step12_prepare.slurm")
printf 'setup_job\t%s\nprep_job\t%s\n' "$SETUP_JOB" "$PREP_JOB"
printf 'After both jobs complete, validate HyPhy help and preparation manifests before submitting arrays.\n'
