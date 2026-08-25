#!/bin/bash

set -euo pipefail

PROJECT=/path/to/mollusc_project
SCRIPT="${PROJECT}/scripts/14_run_class_orthofinder_iqtree.sh"

BIVALVIA_JOB="$(sbatch --parsable --job-name=bivalvia_of_iqtree \
  --export=ALL,GROUP=bivalvia "${SCRIPT}")"
GASTROPODA_JOB="$(sbatch --parsable --job-name=gastropoda_of_iqtree \
  --export=ALL,GROUP=gastropoda "${SCRIPT}")"

echo "Bivalvia OrthoFinder + IQ-TREE: ${BIVALVIA_JOB}"
echo "Gastropoda OrthoFinder + IQ-TREE: ${GASTROPODA_JOB}"
