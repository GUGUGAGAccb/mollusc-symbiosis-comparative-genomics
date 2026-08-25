#!/usr/bin/env bash
#SBATCH --job-name=step12_launch
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step12_launch_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step12_launch_%j.err

set -euo pipefail

PREP_JOB_ID=${1:?Usage: sbatch script PREP_JOB_ID}
ROOT=/path/to/mollusc_project
SCRIPTS="$ROOT/scripts"
RUN="$ROOT/results_or_reports/molecular_evolution/hyphy_busted_absrel_relax_v1/run_${PREP_JOB_ID}"
HYPHY="$ROOT/conda_envs/hyphy_env/bin/hyphy"
JOB_TABLE="$RUN/step12_job_chain.tsv"

test -s "$ROOT/conda_envs/hyphy_env/STEP12_READY"
test -s "$RUN/PREPARED"
test -x "$HYPHY"
for method in busted absrel relax; do
  test -s "$ROOT/conda_envs/hyphy_env/${method}_help.txt"
  grep -q '^alignment' "$ROOT/conda_envs/hyphy_env/${method}_help.txt"
  grep -q '^tree' "$ROOT/conda_envs/hyphy_env/${method}_help.txt"
  grep -q '^output' "$ROOT/conda_envs/hyphy_env/${method}_help.txt"
done
grep -q '^branches' "$ROOT/conda_envs/hyphy_env/busted_help.txt"
grep -q '^branches' "$ROOT/conda_envs/hyphy_env/absrel_help.txt"
grep -q '^test' "$ROOT/conda_envs/hyphy_env/relax_help.txt"
grep -q '^reference' "$ROOT/conda_envs/hyphy_env/relax_help.txt"

for group in bivalvia gastropoda; do
  rows=$(( $(wc -l < "$RUN/manifest_${group}.tsv") - 1 ))
  test "$rows" -eq 30
done

BIVALVIA_JOB=$(sbatch --parsable --array=0-29%3 --job-name=step12_hyphy_bivalvia \
  "$SCRIPTS/29_step12_family_array.slurm" bivalvia "$PREP_JOB_ID")
GASTROPODA_JOB=$(sbatch --parsable --array=0-29%3 --job-name=step12_hyphy_gastropoda \
  "$SCRIPTS/29_step12_family_array.slurm" gastropoda "$PREP_JOB_ID")
SUMMARY_JOB=$(sbatch --parsable --dependency=afterany:"$BIVALVIA_JOB":"$GASTROPODA_JOB" \
  "$SCRIPTS/29_step12_summarise.slurm" "$PREP_JOB_ID")
ARCHIVE_JOB=$(sbatch --parsable --dependency=afterok:"$SUMMARY_JOB" \
  "$SCRIPTS/29_step12_archive.slurm")

printf 'role\tjob_id\ndependency_setup\t%s\ndependency_prepare\t%s\nbivalvia_array\t%s\ngastropoda_array\t%s\nsummary\t%s\narchive\t%s\n' \
  "${SLURM_JOB_DEPENDENCY:-afterok_external}" "$PREP_JOB_ID" "$BIVALVIA_JOB" "$GASTROPODA_JOB" "$SUMMARY_JOB" "$ARCHIVE_JOB" > "$JOB_TABLE"
printf 'Launched Step 12 arrays: bivalvia=%s gastropoda=%s summary=%s archive=%s\n' \
  "$BIVALVIA_JOB" "$GASTROPODA_JOB" "$SUMMARY_JOB" "$ARCHIVE_JOB"
