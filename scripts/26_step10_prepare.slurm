#!/usr/bin/env bash
#SBATCH --job-name=step10_prepare
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step10_prepare_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step10_prepare_%j.err

set -euo pipefail

ROOT=/path/to/mollusc_project
BASE="$ROOT/results_or_reports/candidate_family_phylogenetics/validated_top30_v1"
RUN="$BASE/run_${SLURM_JOB_ID}"
SCRIPT="$ROOT/scripts/26_step10_candidate_gene_phylogenetics.py"
mkdir -p "$RUN" "$ROOT/logs"

test -s "$ROOT/results_or_reports/candidate_sequence_validation/combined_summary/latest/COMPLETE"
python "$SCRIPT" prepare --project-root "$ROOT" --run-dir "$RUN" \
  --top-n 30 --contraction-quota 6 \
  --min-sequences 4 --min-species 3 --min-unique-sequences 3 \
  --per-species-cap 20 --total-cap 180

ln -sfn "run_${SLURM_JOB_ID}" "$BASE/latest_prepared"
printf 'Step 10 inputs prepared: %s\n' "$RUN"
