#!/usr/bin/env bash
#SBATCH --job-name=step08_enrichment
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step08_enrichment_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step08_enrichment_%j.err

set -euo pipefail

ROOT=/path/to/mollusc_project
SCRIPT="$ROOT/scripts/24_step08_functional_enrichment.py"
STEP7="$ROOT/results_or_reports/functional_annotation"
BASE="$ROOT/results_or_reports/functional_enrichment/orthogroup_level_v1"
RUN="$BASE/run_${SLURM_JOB_ID}"
GO_OBO="$ROOT/databases/gene_ontology/go-basic.obo"

mkdir -p "$RUN" "$ROOT/logs"
test -s "$STEP7/bivalvia/eggnog_mapper_v2_1_12/latest/COMPLETE"
test -s "$STEP7/gastropoda/eggnog_mapper_v2_1_12/latest/COMPLETE"
test -s "$GO_OBO"

python "$SCRIPT" \
  --bivalvia "$STEP7/bivalvia/eggnog_mapper_v2_1_12/latest" \
  --gastropoda "$STEP7/gastropoda/eggnog_mapper_v2_1_12/latest" \
  --go-obo "$GO_OBO" \
  --output-dir "$RUN" \
  --min-target-count 3 \
  --min-background-count 5 \
  --min-annotated-target 5 \
  --min-species-set 5

ln -sfn "run_${SLURM_JOB_ID}" "$BASE/latest"
printf 'Step 8 complete: %s\n' "$RUN"
