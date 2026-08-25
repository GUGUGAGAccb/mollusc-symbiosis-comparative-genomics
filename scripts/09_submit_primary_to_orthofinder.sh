#!/bin/bash

set -euo pipefail

SCRATCH=/path/to/mollusc_project
SCRIPTS="${SCRATCH}/scripts"

mkdir -p "${SCRATCH}/logs"

PREP_JOB="$(sbatch --parsable "${SCRIPTS}/04_prepare_primary_proteomes.sh")"
BUSCO_JOB="$(sbatch --parsable --dependency="afterok:${PREP_JOB}" --array=1-14%4 "${SCRIPTS}/05_primary_busco_qc.sh")"
SUMMARY_JOB="$(sbatch --parsable --dependency="afterok:${BUSCO_JOB}" "${SCRIPTS}/06_summarise_primary_busco.sh")"
ORTHOFINDER_JOB="$(sbatch --parsable --dependency="afterok:${SUMMARY_JOB}" "${SCRIPTS}/08_run_orthofinder.sh")"

echo "Primary-proteome preparation: ${PREP_JOB}"
echo "Primary-proteome BUSCO array: ${BUSCO_JOB}"
echo "BUSCO summary: ${SUMMARY_JOB}"
echo "OrthoFinder core analysis: ${ORTHOFINDER_JOB}"
