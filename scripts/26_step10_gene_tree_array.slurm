#!/usr/bin/env bash
#SBATCH --job-name=step10_gene_tree
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step10_gene_tree_%x_%A_%a.out
#SBATCH --error=/path/to/mollusc_project/logs/step10_gene_tree_%x_%A_%a.err

set -euo pipefail

GROUP=${1:?Usage: sbatch --array=0-29 script GROUP}
case "$GROUP" in bivalvia|gastropoda) ;; *) echo "Unsupported group: $GROUP" >&2; exit 2 ;; esac

ROOT=/path/to/mollusc_project
BASE="$ROOT/results_or_reports/candidate_family_phylogenetics/validated_top30_v1"
RUN=$(readlink -f "$BASE/latest_prepared")
MANIFEST="$RUN/selected_${GROUP}.tsv"
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$MANIFEST")
if [[ -z "$LINE" ]]; then
  printf 'No selected family for %s array index %s; task is intentionally empty.\n' "$GROUP" "$SLURM_ARRAY_TASK_ID"
  exit 0
fi

IFS=$'\t' read -r RANK ROW_GROUP ORTHOGROUP INPUT_FASTA SEQUENCE_MAP NSEQ NSPEC REST <<< "$LINE"
test "$ROW_GROUP" = "$GROUP"
test -s "$INPUT_FASTA"
OUT="$RUN/trees/$GROUP/$ORTHOGROUP"
mkdir -p "$OUT"
rm -f "$OUT/FAILED" "$OUT/COMPLETE" "$OUT/TREE_PREFIX" "$OUT/SUPPORT_MODE"

MAFFT=/path/to/mollusc_project/conda_envs/orthofinder_env/bin/mafft
IQTREE=/path/to/mollusc_project/conda_envs/orthofinder_env/bin/iqtree
module purge
module load apps/trimal/1.5

fail_soft() {
  printf '%s\n' "$1" > "$OUT/FAILED"
  printf 'WARNING: %s %s: %s\n' "$GROUP" "$ORTHOGROUP" "$1" >&2
  exit 0
}

if ! "$MAFFT" --auto --thread "${SLURM_CPUS_PER_TASK}" "$INPUT_FASTA" > "$OUT/aligned.faa" 2> "$OUT/mafft.log"; then
  fail_soft "MAFFT failed"
fi
test -s "$OUT/aligned.faa" || fail_soft "MAFFT produced an empty alignment"

if ! trimal -in "$OUT/aligned.faa" -out "$OUT/trimmed.faa" -automated1 > "$OUT/trimal.log" 2>&1; then
  cp "$OUT/aligned.faa" "$OUT/trimmed.faa"
  printf 'trimAl failed; using the untrimmed MAFFT alignment.\n' > "$OUT/TRIM_FALLBACK"
fi

TRIMMED_N=$(grep -c '^>' "$OUT/trimmed.faa" || true)
TRIMMED_LEN=$(awk 'BEGIN{n=0} /^>/{if(n>0){print n; exit} next} {n+=length($0)} END{if(n>0) print n}' "$OUT/trimmed.faa" | head -n 1)
if [[ "$TRIMMED_N" -lt 4 || "${TRIMMED_LEN:-0}" -lt 30 ]]; then
  cp "$OUT/aligned.faa" "$OUT/trimmed.faa"
  printf 'Trimmed alignment was too small; using the untrimmed MAFFT alignment.\n' > "$OUT/TRIM_FALLBACK"
fi

PREFIX_SUPPORTED="$OUT/gene_tree_supported"
MODEL_ARGS=(-m MFP -mset LG,WAG,JTT,Q.INSECT,Q.YEAST -mfreq FU,F -mrate E,I,G,I+G,R -cmax 4)
if "$IQTREE" -s "$OUT/trimmed.faa" "${MODEL_ARGS[@]}" -B 1000 --alrt 1000 -T "${SLURM_CPUS_PER_TASK}" --prefix "$PREFIX_SUPPORTED" --redo > "$OUT/iqtree_supported.stdout" 2> "$OUT/iqtree_supported.stderr"; then
  test -s "$PREFIX_SUPPORTED.treefile" || fail_soft "IQ-TREE support run returned without a treefile"
  printf '%s\n' "$PREFIX_SUPPORTED" > "$OUT/TREE_PREFIX"
  printf 'ufboot1000_alrt1000\n' > "$OUT/SUPPORT_MODE"
else
  PREFIX_FALLBACK="$OUT/gene_tree_ml_only"
  printf 'Supported IQ-TREE run failed; attempting ML-only fallback.\n' > "$OUT/IQTREE_SUPPORT_FALLBACK"
  if ! "$IQTREE" -s "$OUT/trimmed.faa" "${MODEL_ARGS[@]}" -T "${SLURM_CPUS_PER_TASK}" --prefix "$PREFIX_FALLBACK" --redo > "$OUT/iqtree_fallback.stdout" 2> "$OUT/iqtree_fallback.stderr"; then
    fail_soft "Both supported and ML-only IQ-TREE runs failed"
  fi
  test -s "$PREFIX_FALLBACK.treefile" || fail_soft "ML-only IQ-TREE fallback returned without a treefile"
  printf '%s\n' "$PREFIX_FALLBACK" > "$OUT/TREE_PREFIX"
  printf 'ml_only_fallback\n' > "$OUT/SUPPORT_MODE"
fi

printf 'group\torthogroup\trank\tinput_sequences\tinput_species\n%s\t%s\t%s\t%s\t%s\n' \
  "$GROUP" "$ORTHOGROUP" "$RANK" "$NSEQ" "$NSPEC" > "$OUT/task_metadata.tsv"
printf 'Step 10 gene tree complete\n' > "$OUT/COMPLETE"
printf 'Completed %s %s (%s sequences, %s species)\n' "$GROUP" "$ORTHOGROUP" "$NSEQ" "$NSPEC"
