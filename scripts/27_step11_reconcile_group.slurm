#!/usr/bin/env bash
#SBATCH --job-name=step11_reconcile
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/step11_reconcile_%x_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/step11_reconcile_%x_%j.err

set -euo pipefail

GROUP=${1:?Usage: sbatch script GROUP}
case "$GROUP" in bivalvia|gastropoda) ;; *) echo "Unsupported group: $GROUP" >&2; exit 2 ;; esac

ROOT=/path/to/mollusc_project
STEP10="$ROOT/results_or_reports/candidate_family_phylogenetics/validated_top30_v1/latest"
BASE="$ROOT/results_or_reports/gene_tree_species_tree_reconciliation/$GROUP/dl_lca_v1"
RUN="$BASE/run_${SLURM_JOB_ID}"
SCRIPT="$ROOT/scripts/27_step11_reconciliation.py"
PYTHON="$ROOT/conda_envs/orthofinder_env/bin/python"

if [[ "$GROUP" == bivalvia ]]; then
  CLASSROOT="$ROOT/orthofinder_results/by_class/bivalvia_v3_1_5/Results_Aug04"
  SPECIES_TREE="$CLASSROOT/Species_Tree/SpeciesTree_rooted.txt"
  IQTREE="$CLASSROOT/Species_Tree_IQTree/SpeciesTree_IQTREE.treefile"
else
  CLASSROOT="$ROOT/orthofinder_results/by_class/gastropoda_v3_1_5/Results_Aug04"
  SPECIES_TREE="$CLASSROOT/Species_Tree/SpeciesTree_rooted.txt"
  IQTREE="$CLASSROOT/Species_Tree_IQTree_restricted_mfp_v1/SpeciesTree_IQTREE.treefile"
fi

mkdir -p "$RUN" "$ROOT/logs"
test -s "$STEP10/COMPLETE"
test -s "$SPECIES_TREE"
test -s "$IQTREE"

"$PYTHON" "$SCRIPT" reconcile-group \
  --group "$GROUP" \
  --step10 "$STEP10" \
  --species-tree "$SPECIES_TREE" \
  --iqtree-species-tree "$IQTREE" \
  --output-dir "$RUN" \
  --duplication-cost 1 --loss-cost 1

ln -sfn "run_${SLURM_JOB_ID}" "$BASE/latest"
printf 'Step 11 %s reconciliation complete: %s\n' "$GROUP" "$RUN"
