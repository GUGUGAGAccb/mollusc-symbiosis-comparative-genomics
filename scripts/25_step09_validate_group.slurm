#!/usr/bin/env bash
#SBATCH --job-name=step09_validate
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step09_validate_%x_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step09_validate_%x_%j.err

set -euo pipefail

GROUP=${1:?Usage: sbatch script GROUP}
case "$GROUP" in bivalvia|gastropoda) ;; *) echo "Unsupported group: $GROUP" >&2; exit 2 ;; esac

ROOT=/path/to/mollusc_project
STEP7="$ROOT/results_or_reports/functional_annotation/$GROUP/eggnog_mapper_v2_1_12/latest"
PFAM="$ROOT/databases/pfam/current"
BASE="$ROOT/results_or_reports/candidate_sequence_validation/$GROUP/pfam_hmmer_v1"
RUN="$BASE/run_${SLURM_JOB_ID}"
SCRIPT="$ROOT/scripts/25_step09_candidate_sequence_validation.py"
TAXDB=${EGGNOG_TAXDB:-/path/to/eggnog.taxa.db}
mkdir -p "$RUN" "$ROOT/logs"

test -s "$STEP7/COMPLETE"
test -s "$PFAM/READY"
test -s "$TAXDB"
module purge
module load apps/hmmer/3.4

hmmscan --cpu "${SLURM_CPUS_PER_TASK}" --cut_ga --noali \
  --domtblout "$RUN/direct_pfam.domtblout" \
  "$PFAM/Pfam-A.hmm" "$STEP7/annotation_input_unique.faa" \
  > "$RUN/hmmscan.out"

test -s "$RUN/direct_pfam.domtblout"
python "$SCRIPT" validate \
  --group "$GROUP" \
  --fasta "$STEP7/annotation_input_unique.faa" \
  --annotations "$STEP7/protein_functional_annotations.tsv" \
  --family-table "$STEP7/high_confidence_family_functional_annotations.tsv" \
  --domtblout "$RUN/direct_pfam.domtblout" \
  --taxonomy-db "$TAXDB" \
  --output-dir "$RUN"

ln -sfn "run_${SLURM_JOB_ID}" "$BASE/latest"
printf 'Step 9 %s validation complete: %s\n' "$GROUP" "$RUN"
