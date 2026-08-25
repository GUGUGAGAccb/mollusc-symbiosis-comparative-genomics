#!/bin/bash
#SBATCH --job-name=cafe5_group
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step05_cafe5_%x_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step05_cafe5_%x_%j.err

set -euo pipefail

: "${GROUP:?Submit with --export=ALL,GROUP=bivalvia or GROUP=gastropoda}"
if [[ "$GROUP" != "bivalvia" && "$GROUP" != "gastropoda" ]]; then
    echo "Unsupported GROUP=$GROUP" >&2
    exit 2
fi

PROJECT=/path/to/mollusc_project
ORTHO_ENV="$PROJECT/conda_envs/orthofinder_env"
CAFE_ENV="$PROJECT/conda_envs/cafe5_env"
PREP_SCRIPT="$PROJECT/scripts/17_step05_prepare_cafe.py"
BASE="$PROJECT/orthofinder_results/by_class/${GROUP}_v3_1_5/Results_Aug04"
GENE_COUNT="$BASE/Orthogroups/Orthogroups.GeneCount.tsv"
if [[ "$GROUP" == "gastropoda" ]]; then
    IQTREE_PREFIX="$BASE/Species_Tree_IQTree_restricted_mfp_v1/SpeciesTree_IQTREE"
else
    IQTREE_PREFIX="$BASE/Species_Tree_IQTree/SpeciesTree_IQTREE"
fi
REFERENCE_TREE="$BASE/Species_Tree/SpeciesTree_rooted.txt"
OUTPUT_ROOT="$PROJECT/results_or_reports/gene_family_evolution/$GROUP/relative_time_screen_v1"
RUN_DIR="$OUTPUT_ROOT/run_${SLURM_JOB_ID}"

mkdir -p "$OUTPUT_ROOT" "$PROJECT/logs"

for required in "$GENE_COUNT" "$IQTREE_PREFIX.treefile" "$IQTREE_PREFIX.contree" "$IQTREE_PREFIX.iqtree" "$REFERENCE_TREE"; do
    if [[ ! -s "$required" ]]; then
        echo "Required input missing or empty: $required" >&2
        exit 3
    fi
done
test -x "$ORTHO_ENV/bin/python"
test -x "$CAFE_ENV/bin/cafe5"

"$ORTHO_ENV/bin/python" "$PREP_SCRIPT" \
    --group "$GROUP" \
    --gene-count "$GENE_COUNT" \
    --iqtree-tree "$IQTREE_PREFIX.treefile" \
    --reference-rooted-tree "$REFERENCE_TREE" \
    --root-age 100 \
    --max-family-size 99 \
    --output "$RUN_DIR"

TREE="$RUN_DIR/SpeciesTree_IQTREE.rooted_relative_ultrametric.tre"
FAMILIES="$RUN_DIR/cafe_input_primary.tsv"

"$CAFE_ENV/bin/cafe5" \
    --cores "$SLURM_CPUS_PER_TASK" \
    --tree "$TREE" \
    --infile "$FAMILIES" \
    --poisson \
    --error_model \
    --pvalue 0.05 \
    --output_prefix "$RUN_DIR/cafe_base_error"

ERROR_MODEL=$(find "$RUN_DIR/cafe_base_error" -type f -iname '*error*model*' -print -quit)
if [[ -z "$ERROR_MODEL" || ! -s "$ERROR_MODEL" ]]; then
    echo "CAFE5 did not produce an error-model file" >&2
    exit 4
fi

"$CAFE_ENV/bin/cafe5" \
    --cores "$SLURM_CPUS_PER_TASK" \
    --tree "$TREE" \
    --infile "$FAMILIES" \
    --poisson \
    --n_gamma_cats 2 \
    -e"$ERROR_MODEL" \
    --pvalue 0.05 \
    --output_prefix "$RUN_DIR/cafe_gamma_k2"

find "$RUN_DIR/cafe_base_error" "$RUN_DIR/cafe_gamma_k2" -type f -printf '%p\t%s\n' | sort > "$RUN_DIR/cafe_output_manifest.tsv"
date --iso-8601=seconds > "$RUN_DIR/COMPLETE"
ln -sfn "$RUN_DIR" "$OUTPUT_ROOT/latest"

echo "STEP05_COMPLETE=$RUN_DIR"
