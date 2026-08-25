#!/bin/bash
#SBATCH --job-name=validate_cafe_prep
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/validate_cafe_prep_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/validate_cafe_prep_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
BASE="$PROJECT/orthofinder_results/by_class/bivalvia_v3_1_5/Results_Aug04"
OUT="$PROJECT/tmp/step05_prepare_validation_${SLURM_JOB_ID}"

mkdir -p "$PROJECT/tmp"
"$PROJECT/conda_envs/orthofinder_env/bin/python" \
    "$PROJECT/scripts/17_step05_prepare_cafe.py" \
    --group bivalvia \
    --gene-count "$BASE/Orthogroups/Orthogroups.GeneCount.tsv" \
    --iqtree-tree "$BASE/Species_Tree_IQTree/SpeciesTree_IQTREE.treefile" \
    --reference-rooted-tree "$BASE/Species_Tree/SpeciesTree_rooted.txt" \
    --root-age 100 \
    --max-family-size 99 \
    --output "$OUT"

test -s "$OUT/SpeciesTree_IQTREE.rooted_relative_ultrametric.tre"
test -s "$OUT/cafe_input_primary.tsv"
test -s "$OUT/cafe_input_summary.tsv"
cat "$OUT/cafe_input_summary.tsv"
echo "VALIDATION_DIR=$OUT"
