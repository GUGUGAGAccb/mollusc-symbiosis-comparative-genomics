#!/bin/bash

set -euo pipefail

PROJECT=/path/to/mollusc_project
DOWNLOAD_JOB="$(sbatch --parsable "${PROJECT}/scripts/10_download_fragum_sueziense.sh")"
BUSCO_JOB="$(sbatch --parsable --dependency="afterok:${DOWNLOAD_JOB}" \
  "${PROJECT}/scripts/11_busco_fragum_sueziense.sh")"

echo "Fragum_sueziense download/preparation: ${DOWNLOAD_JOB}"
echo "Fragum_sueziense BUSCO: ${BUSCO_JOB}"
