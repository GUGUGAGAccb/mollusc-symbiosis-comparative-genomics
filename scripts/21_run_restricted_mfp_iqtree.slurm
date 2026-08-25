#!/bin/bash
#SBATCH --job-name=iqtree_restricted
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=3-00:00:00
#SBATCH --mem=16G
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/iqtree_restricted_%x_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/iqtree_restricted_%x_%j.err

set -euo pipefail

: "${ANALYSIS:?Submit with --export=ALL,ANALYSIS=combined or ANALYSIS=gastropoda}"

PROJECT=/path/to/mollusc_project
CONDA_BASE=/path/to/miniconda3
CONDA_ENV="${PROJECT}/conda_envs/orthofinder_env"

case "${ANALYSIS}" in
  combined)
    BASE="${PROJECT}/orthofinder_results/core_complete80_v3_1_5/Results_Aug04"
    EXPECTED_SPECIES=14
    ;;
  gastropoda)
    BASE="${PROJECT}/orthofinder_results/by_class/gastropoda_v3_1_5/Results_Aug04"
    EXPECTED_SPECIES=8
    ;;
  *)
    echo "Unsupported ANALYSIS=${ANALYSIS}" >&2
    exit 2
    ;;
esac

INPUT="${BASE}/Species_Tree_IQTree/SpeciesTreeAlignment.species_names.fa"
OUTPUT_DIR="${BASE}/Species_Tree_IQTree_restricted_mfp_v1"
PREFIX="${OUTPUT_DIR}/SpeciesTree_IQTREE"

[[ -s "${INPUT}" ]] || { echo "Missing named alignment: ${INPUT}" >&2; exit 3; }
N_SEQUENCES="$(grep -c '^>' "${INPUT}")"
N_UNIQUE="$(grep '^>' "${INPUT}" | sort -u | wc -l)"
[[ "${N_SEQUENCES}" -eq "${EXPECTED_SPECIES}" && "${N_UNIQUE}" -eq "${EXPECTED_SPECIES}" ]] || {
  echo "Expected ${EXPECTED_SPECIES} unique species, found ${N_SEQUENCES}/${N_UNIQUE}" >&2
  exit 4
}

mkdir -p "${OUTPUT_DIR}" "${PROJECT}/logs"
if [[ -s "${PREFIX}.treefile" || -s "${PREFIX}.iqtree" ]]; then
  echo "Refusing to overwrite a completed restricted IQ-TREE run: ${PREFIX}" >&2
  exit 5
fi

printf '%s\n' \
  "analysis=${ANALYSIS}" \
  "input=${INPUT}" \
  "model_search=MFP" \
  "model_set=Q.INSECT,Q.YEAST,JTT,LG,WAG" \
  "frequencies=FU,F" \
  "rate_models=E,I,G,I+G,R" \
  "max_rate_categories=6" \
  "ufboot=1000" \
  "shalrt=1000" \
  "threads=${SLURM_CPUS_PER_TASK}" \
  "seed=12345" > "${OUTPUT_DIR}/RUN_PARAMETERS.txt"

export CONDA_ENVS_PATH="${PROJECT}/conda_envs"
export CONDA_PKGS_DIRS="${PROJECT}/conda_pkgs"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

iqtree3 \
  -s "${INPUT}" \
  -st AA \
  -m MFP \
  -mset Q.INSECT,Q.YEAST,JTT,LG,WAG \
  -mfreq FU,F \
  -mrate E,I,G,I+G,R \
  -cmax 6 \
  -B 1000 \
  --alrt 1000 \
  -T "${SLURM_CPUS_PER_TASK}" \
  -seed 12345 \
  --prefix "${PREFIX}"

for required in "${PREFIX}.treefile" "${PREFIX}.iqtree" "${PREFIX}.contree"; do
  [[ -s "${required}" ]] || { echo "Missing final IQ-TREE output: ${required}" >&2; exit 6; }
done

date --iso-8601=seconds > "${OUTPUT_DIR}/COMPLETE"
echo "RESTRICTED_IQTREE_COMPLETE=${OUTPUT_DIR}"
