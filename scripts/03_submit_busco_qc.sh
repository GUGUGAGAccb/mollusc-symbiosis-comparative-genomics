#!/bin/bash

set -euo pipefail

SCRATCH=/path/to/mollusc_project
PROTEOME_DIR="${SCRATCH}/orthofinder_proteomes"
SCRIPT_DIR="${SCRATCH}/scripts"

mkdir -p "${SCRATCH}/logs"

N_PROTEOMES="$(
  find "${PROTEOME_DIR}" -maxdepth 1 -type f -name '*.faa' -size +0c | wc -l
)"

if (( N_PROTEOMES == 0 )); then
  echo "ERROR: no non-empty .faa files found in ${PROTEOME_DIR}" >&2
  exit 1
fi

QC_JOB_ID="$(
  sbatch --parsable --array="1-${N_PROTEOMES}%4" "${SCRIPT_DIR}/01_proteome_busco_qc.sh"
)"

SUMMARY_JOB_ID="$(
  sbatch --parsable --dependency="afterok:${QC_JOB_ID}" "${SCRIPT_DIR}/02_summarise_busco.sh"
)"

echo "Submitted BUSCO array job: ${QC_JOB_ID} (1-${N_PROTEOMES}%4)"
echo "Submitted dependent summary job: ${SUMMARY_JOB_ID}"
