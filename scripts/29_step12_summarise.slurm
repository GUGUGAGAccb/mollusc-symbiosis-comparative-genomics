#!/usr/bin/env bash
#SBATCH --job-name=step12_summarise
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step12_summarise_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step12_summarise_%j.err

set -euo pipefail

PREP_JOB_ID=${1:?Usage: sbatch script PREP_JOB_ID}
ROOT=/path/to/mollusc_project
BASE="$ROOT/results_or_reports/molecular_evolution/hyphy_busted_absrel_relax_v1"
RUN="$BASE/run_${PREP_JOB_ID}"
PYTHON="$ROOT/conda_envs/orthofinder_env/bin/python"
SCRIPT="$ROOT/scripts/29_step12_molecular_evolution.py"

test -s "$RUN/PREPARED"
"$PYTHON" "$SCRIPT" summarise --run-dir "$RUN"
test -s "$RUN/COMPLETE"
ln -sfn "run_${PREP_JOB_ID}" "$BASE/latest"
printf 'Step 12 summary complete: %s\n' "$RUN"

