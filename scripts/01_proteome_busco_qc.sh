#!/bin/bash
#SBATCH --job-name=proteome_busco
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --array=1-14%4
#SBATCH --output=/path/to/mollusc_project/logs/busco_%A_%a.out
#SBATCH --error=/path/to/mollusc_project/logs/busco_%A_%a.err

set -euo pipefail

SCRATCH=/path/to/mollusc_project
PROTEOME_DIR="${SCRATCH}/orthofinder_proteomes"
BUSCO_DB_DIR="${SCRATCH}/busco_downloads"
BUSCO_OUT_DIR="${SCRATCH}/busco_proteomes"
QC_DIR="${SCRATCH}/proteome_qc"
CONDA_BASE=/path/to/miniconda3
CONDA_ENV="${SCRATCH}/conda_envs/busco_env"

export CONDA_ENVS_PATH="${SCRATCH}/conda_envs"
export CONDA_PKGS_DIRS="${SCRATCH}/conda_pkgs"
export LC_ALL=C

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p \
  "${BUSCO_OUT_DIR}" \
  "${QC_DIR}/seqkit_stats" \
  "${QC_DIR}/duplicate_ids" \
  "${QC_DIR}/invalid_characters" \
  "${QC_DIR}/busco_summaries" \
  "${SCRATCH}/logs"

mapfile -t PROTEOMES < <(
  find "${PROTEOME_DIR}" -maxdepth 1 -type f -name '*.faa' -size +0c | sort
)

N_PROTEOMES="${#PROTEOMES[@]}"
if (( N_PROTEOMES == 0 )); then
  echo "ERROR: no non-empty .faa files found in ${PROTEOME_DIR}" >&2
  exit 1
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
if (( TASK_ID < 1 || TASK_ID > N_PROTEOMES )); then
  echo "ERROR: array task ${TASK_ID} is outside 1-${N_PROTEOMES}" >&2
  exit 1
fi

PROTEOME="${PROTEOMES[$((TASK_ID - 1))]}"
BASENAME="$(basename "${PROTEOME}")"
SPECIES="${BASENAME%.faa}"

echo "Species: ${SPECIES}"
echo "Proteome: ${PROTEOME}"
echo "Threads: ${THREADS}"

if [[ "$(grep -m1 -v '^[[:space:]]*$' "${PROTEOME}" | cut -c1)" != '>' ]]; then
  echo "ERROR: file does not appear to be FASTA: ${PROTEOME}" >&2
  exit 1
fi

N_PROTEINS="$(grep -c '^>' "${PROTEOME}")"
if (( N_PROTEINS == 0 )); then
  echo "ERROR: no protein sequences found in ${PROTEOME}" >&2
  exit 1
fi

seqkit stats --tabular --all --threads "${THREADS}" "${PROTEOME}" \
  > "${QC_DIR}/seqkit_stats/${SPECIES}.seqkit_stats.tsv"

DUPLICATE_FILE="${QC_DIR}/duplicate_ids/${SPECIES}.duplicate_ids.txt"
grep '^>' "${PROTEOME}" | sed 's/^>//' | awk '{print $1}' | sort | uniq -d \
  > "${DUPLICATE_FILE}"

INVALID_FILE="${QC_DIR}/invalid_characters/${SPECIES}.invalid_characters.txt"
grep -v '^>' "${PROTEOME}" \
  | tr -d '\n\r[:space:]' \
  | fold -w1 \
  | tr '[:lower:]' '[:upper:]' \
  | grep -v -E '^[ABCDEFGHIKLMNPQRSTVWYXBZJUO*]$' \
  | sort \
  | uniq -c \
  > "${INVALID_FILE}" || true

SPECIES_OUT_DIR="${BUSCO_OUT_DIR}/${SPECIES}"
mkdir -p "${SPECIES_OUT_DIR}"

busco \
  --in "${PROTEOME}" \
  --out "${SPECIES}_metazoa_odb10" \
  --out_path "${SPECIES_OUT_DIR}" \
  --mode proteins \
  --lineage_dataset "${BUSCO_DB_DIR}/lineages/metazoa_odb10" \
  --cpu "${THREADS}" \
  --offline \
  --force

SUMMARY_SOURCE="$(
  find "${SPECIES_OUT_DIR}/${SPECIES}_metazoa_odb10" \
    -maxdepth 1 -type f -name 'short_summary*.txt' | head -n 1
)"

if [[ -z "${SUMMARY_SOURCE}" ]]; then
  echo "ERROR: BUSCO short summary not found for ${SPECIES}" >&2
  exit 1
fi

cp "${SUMMARY_SOURCE}" "${QC_DIR}/busco_summaries/${SPECIES}.busco_short_summary.txt"

echo "Completed ${SPECIES}: ${N_PROTEINS} proteins"
