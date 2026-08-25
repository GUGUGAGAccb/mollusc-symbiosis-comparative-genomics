#!/bin/bash
#SBATCH --job-name=iqtree_species
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=3-00:00:00
#SBATCH --mem=80G
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/iqtree_species_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/iqtree_species_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
OF_RUN="${PROJECT}/orthofinder_results/core_complete80_v3_1_5/Results_Aug04"
INPUT_ALIGNMENT="${OF_RUN}/WorkingDirectory/Alignments_ids/SpeciesTreeAlignment.fa"
SPECIES_IDS="${OF_RUN}/WorkingDirectory/SpeciesIDs.txt"
OUTPUT_DIR="${OF_RUN}/Species_Tree_IQTree"
RENAMED_ALIGNMENT="${OUTPUT_DIR}/SpeciesTreeAlignment.species_names.fa"
PREFIX="${OUTPUT_DIR}/SpeciesTree_IQTREE"
CONDA_BASE=/path/to/miniconda3
CONDA_ENV="${PROJECT}/conda_envs/orthofinder_env"

[[ -s "${INPUT_ALIGNMENT}" ]] || { echo "Missing alignment: ${INPUT_ALIGNMENT}" >&2; exit 1; }
[[ -s "${SPECIES_IDS}" ]] || { echo "Missing species mapping: ${SPECIES_IDS}" >&2; exit 1; }

mkdir -p "${OUTPUT_DIR}" "${PROJECT}/logs"
if compgen -G "${PREFIX}.*" > /dev/null; then
  echo "Refusing to overwrite an existing IQ-TREE run: ${PREFIX}.*" >&2
  exit 1
fi

awk '
  FNR == NR {
    split($0, fields, ": ")
    species = fields[2]
    sub(/\.faa$/, "", species)
    names[fields[1]] = species
    next
  }
  /^>/ {
    id = substr($0, 2)
    if (!(id in names)) {
      print "Unmapped species ID: " id > "/dev/stderr"
      exit 2
    }
    print ">" names[id]
    next
  }
  { print }
' "${SPECIES_IDS}" "${INPUT_ALIGNMENT}" > "${RENAMED_ALIGNMENT}"

N_SEQUENCES="$(grep -c '^>' "${RENAMED_ALIGNMENT}")"
N_UNIQUE="$(grep '^>' "${RENAMED_ALIGNMENT}" | sort -u | wc -l)"
[[ "${N_SEQUENCES}" -eq 14 && "${N_UNIQUE}" -eq 14 ]] || {
  echo "Expected 14 unique species, found ${N_SEQUENCES} sequences and ${N_UNIQUE} unique names" >&2
  exit 1
}

export CONDA_ENVS_PATH="${PROJECT}/conda_envs"
export CONDA_PKGS_DIRS="${PROJECT}/conda_pkgs"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

iqtree3 \
  -s "${RENAMED_ALIGNMENT}" \
  -st AA \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T "${SLURM_CPUS_PER_TASK}" \
  -seed 12345 \
  --prefix "${PREFIX}"

[[ -s "${PREFIX}.treefile" ]] || { echo "IQ-TREE treefile was not produced" >&2; exit 1; }
[[ -s "${PREFIX}.iqtree" ]] || { echo "IQ-TREE report was not produced" >&2; exit 1; }

echo "IQ-TREE species tree completed"
echo "Alignment: ${RENAMED_ALIGNMENT}"
echo "Tree: ${PREFIX}.treefile"
echo "Report: ${PREFIX}.iqtree"
