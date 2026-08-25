#!/bin/bash
#SBATCH --job-name=class_of_iqtree
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=3-00:00:00
#SBATCH --mem=16G
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/class_of_iqtree_%x_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/class_of_iqtree_%x_%j.err

set -euo pipefail

: "${GROUP:?Submit with --export=ALL,GROUP=bivalvia or GROUP=gastropoda}"

PROJECT=/path/to/mollusc_project
SOURCE_DIR="${PROJECT}/orthofinder_primary_proteomes"
INPUT_PARENT="${PROJECT}/orthofinder_by_class_inputs"
RESULTS_PARENT="${PROJECT}/orthofinder_results/by_class"
INPUT_DIR="${INPUT_PARENT}/${GROUP}"
RESULTS_DIR="${RESULTS_PARENT}/${GROUP}_v3_1_5"
CONDA_BASE=/path/to/miniconda3
CONDA_ENV="${PROJECT}/conda_envs/orthofinder_env"

case "${GROUP}" in
  bivalvia)
    SPECIES=(
      Acanthocardia_echinata
      Americardia_media
      Cerastoderma_edule
      Fragum_sueziense
      Tridacna_gigas
      Tridacna_maxima
    )
    ;;
  gastropoda)
    SPECIES=(
      Aplysia_californica
      Berghia_stephanieae
      Bullacta_exarata
      Elysia_chlorotica
      Elysia_crispata
      Elysia_marginata
      Onchidella_celtica
      Plakobranchus_ocellatus
    )
    ;;
  *)
    echo "Unsupported GROUP: ${GROUP}" >&2
    exit 2
    ;;
esac

if [[ -e "${RESULTS_DIR}" ]]; then
  echo "Refusing to overwrite existing results: ${RESULTS_DIR}" >&2
  exit 1
fi

mkdir -p "${INPUT_PARENT}" "${RESULTS_PARENT}" "${PROJECT}/logs"
if [[ -e "${INPUT_DIR}" ]]; then
  echo "Refusing to reuse existing input directory: ${INPUT_DIR}" >&2
  exit 1
fi
mkdir "${INPUT_DIR}"

for species in "${SPECIES[@]}"; do
  source_fasta="${SOURCE_DIR}/${species}.faa"
  [[ -s "${source_fasta}" ]] || { echo "Missing proteome: ${source_fasta}" >&2; exit 1; }
  ln -s "${source_fasta}" "${INPUT_DIR}/${species}.faa"
done

EXPECTED="${#SPECIES[@]}"
FOUND="$(find "${INPUT_DIR}" -maxdepth 1 -type l -name '*.faa' | wc -l)"
[[ "${FOUND}" -eq "${EXPECTED}" ]] || {
  echo "Expected ${EXPECTED} proteomes, found ${FOUND}" >&2
  exit 1
}

export CONDA_ENVS_PATH="${PROJECT}/conda_envs"
export CONDA_PKGS_DIRS="${PROJECT}/conda_pkgs"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
export TMPDIR="${PROJECT}/tmp/${GROUP}_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

echo "Group: ${GROUP}"
echo "Species count: ${EXPECTED}"
printf '%s\n' "${SPECIES[@]}"
orthofinder --version

orthofinder \
  -f "${INPUT_DIR}" \
  -t "${SLURM_CPUS_PER_TASK}" \
  -a "${SLURM_CPUS_PER_TASK}" \
  -M msa \
  -S diamond \
  -A famsa \
  -T fasttree \
  -o "${RESULTS_DIR}"

mapfile -t RUN_DIRS < <(find "${RESULTS_DIR}" -mindepth 1 -maxdepth 1 -type d -name 'Results_*' | sort)
[[ "${#RUN_DIRS[@]}" -eq 1 ]] || {
  echo "Expected one OrthoFinder Results_* directory, found ${#RUN_DIRS[@]}" >&2
  exit 1
}
RUN_DIR="${RUN_DIRS[0]}"
INPUT_ALIGNMENT="${RUN_DIR}/WorkingDirectory/Alignments_ids/SpeciesTreeAlignment.fa"
SPECIES_IDS="${RUN_DIR}/WorkingDirectory/SpeciesIDs.txt"
IQTREE_DIR="${RUN_DIR}/Species_Tree_IQTree"
RENAMED_ALIGNMENT="${IQTREE_DIR}/SpeciesTreeAlignment.species_names.fa"
PREFIX="${IQTREE_DIR}/SpeciesTree_IQTREE"

[[ -s "${INPUT_ALIGNMENT}" ]] || { echo "Missing alignment: ${INPUT_ALIGNMENT}" >&2; exit 1; }
[[ -s "${SPECIES_IDS}" ]] || { echo "Missing species mapping: ${SPECIES_IDS}" >&2; exit 1; }
mkdir "${IQTREE_DIR}"

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
[[ "${N_SEQUENCES}" -eq "${EXPECTED}" && "${N_UNIQUE}" -eq "${EXPECTED}" ]] || {
  echo "Expected ${EXPECTED} unique species in alignment, found ${N_SEQUENCES}/${N_UNIQUE}" >&2
  exit 1
}

iqtree3 \
  -s "${RENAMED_ALIGNMENT}" \
  -st AA \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T "${SLURM_CPUS_PER_TASK}" \
  -seed 12345 \
  --prefix "${PREFIX}"

[[ -s "${PREFIX}.treefile" ]] || { echo "Missing IQ-TREE treefile" >&2; exit 1; }
[[ -s "${PREFIX}.iqtree" ]] || { echo "Missing IQ-TREE report" >&2; exit 1; }

ln -s "$(basename "${RUN_DIR}")" "${RESULTS_DIR}/latest"

echo "Completed group: ${GROUP}"
echo "OrthoFinder results: ${RUN_DIR}"
echo "IQ-TREE tree: ${PREFIX}.treefile"

