#!/bin/bash
#SBATCH --job-name=primary_proteomes
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/primary_proteomes_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/primary_proteomes_%j.err

set -euo pipefail

SCRATCH=/path/to/mollusc_project
SCRIPT_DIR="${SCRATCH}/scripts"
INPUT_DIR="${SCRATCH}/orthofinder_proteomes"
OUTPUT_DIR="${SCRATCH}/orthofinder_primary_proteomes"
MAPPING_DIR="${SCRATCH}/primary_isoform_mappings"
SUMMARY="${SCRATCH}/primary_proteome_selection.tsv"
PYTHON="${SCRATCH}/conda_envs/orthofinder_env/bin/python3"

mkdir -p "${OUTPUT_DIR}" "${MAPPING_DIR}" "${SCRATCH}/logs"

"${PYTHON}" "${SCRIPT_DIR}/04_select_primary_isoforms.py" \
  --input-dir "${INPUT_DIR}" \
  --output-dir "${OUTPUT_DIR}" \
  --mapping-dir "${MAPPING_DIR}" \
  --summary "${SUMMARY}"

for FASTA in "${OUTPUT_DIR}"/*.faa; do
  [[ -s "${FASTA}" ]] || { echo "ERROR: empty output ${FASTA}" >&2; exit 1; }
  IDS="$(grep -c '^>' "${FASTA}")"
  UNIQUE_IDS="$(grep '^>' "${FASTA}" | sed 's/^>//' | awk '{print $1}' | sort -u | wc -l)"
  [[ "${IDS}" -eq "${UNIQUE_IDS}" ]] || { echo "ERROR: duplicate IDs in ${FASTA}" >&2; exit 1; }
done

column -t -s $'\t' "${SUMMARY}"
