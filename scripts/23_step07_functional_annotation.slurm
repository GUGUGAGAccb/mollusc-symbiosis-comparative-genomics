#!/usr/bin/env bash
#SBATCH --job-name=step07_eggnog
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step07_eggnog_%x_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step07_eggnog_%x_%j.err

set -euo pipefail

GROUP=${1:?Usage: sbatch script GROUP}
case "$GROUP" in
  bivalvia|gastropoda) ;;
  *) echo "Unsupported group: $GROUP" >&2; exit 2 ;;
esac

ROOT=/path/to/mollusc_project
STEP6="$ROOT/results_or_reports/candidate_gene_families/base_model_v1/latest/$GROUP"
BASE="$ROOT/results_or_reports/functional_annotation/$GROUP/eggnog_mapper_v2_1_12"
RUN="$BASE/run_${SLURM_JOB_ID}"
RAW="$RUN/raw"
WORKDIR="${TMPDIR:-/tmp}/step07_${GROUP}_${SLURM_JOB_ID}"
SCRIPT="$ROOT/scripts/23_step07_functional_annotation.py"

mkdir -p "$RAW" "$WORKDIR" "$ROOT/logs"
trap 'rm -rf "$WORKDIR"' EXIT

module purge
module load apps/eggnog-mapper/2.1.12

python "$SCRIPT" prepare \
  --group "$GROUP" \
  --high-confidence-members "$STEP6/high_confidence_family_all_members.faa" \
  --candidate-representatives "$STEP6/annotation_representatives_longest_per_species.faa" \
  --output-fasta "$RUN/annotation_input_unique.faa" \
  --inventory "$RUN/annotation_query_inventory.tsv"

emapper.py \
  -i "$RUN/annotation_input_unique.faa" \
  --itype proteins \
  --output "${GROUP}_candidate_families" \
  --output_dir "$RAW" \
  --temp_dir "$WORKDIR" \
  --cpu "${SLURM_CPUS_PER_TASK}" \
  --tax_scope eukaryota \
  --tax_scope_mode inner_narrowest \
  --target_orthologs all \
  --go_evidence non-electronic \
  --pfam_realign none \
  --sensmode sensitive \
  --dmnd_iterate yes \
  --dmnd_algo auto \
  --evalue 0.001 \
  --seed_ortholog_evalue 0.001 \
  --override

ANNOT="$RAW/${GROUP}_candidate_families.emapper.annotations"
SEEDS="$RAW/${GROUP}_candidate_families.emapper.seed_orthologs"
test -s "$ANNOT"

python "$SCRIPT" summarise \
  --group "$GROUP" \
  --inventory "$RUN/annotation_query_inventory.tsv" \
  --annotations "$ANNOT" \
  --seed-orthologs "$SEEDS" \
  --ranked-families "$STEP6/candidate_families_ranked.tsv" \
  --terminal-events "$STEP6/terminal_candidate_events.tsv" \
  --output-dir "$RUN"

ln -sfn "run_${SLURM_JOB_ID}" "$BASE/latest"
printf 'Step 7 %s functional annotation complete: %s\n' "$GROUP" "$RUN"
