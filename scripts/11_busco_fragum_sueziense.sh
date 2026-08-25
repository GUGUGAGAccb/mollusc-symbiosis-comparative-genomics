#!/bin/bash
#SBATCH --job-name=busco_fsue
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/busco_fragum_sueziense_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/busco_fragum_sueziense_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
SPECIES=Fragum_sueziense
PROTEOME="${PROJECT}/orthofinder_primary_proteomes/${SPECIES}.faa"
BUSCO_DB="${PROJECT}/busco_downloads/lineages/metazoa_odb10"
BUSCO_OUT="${PROJECT}/busco_primary_proteomes/${SPECIES}"
CONDA_BASE=/path/to/miniconda3
CONDA_ENV="${PROJECT}/conda_envs/busco_env"

[[ -s "${PROTEOME}" ]] || { echo "Missing representative proteome: ${PROTEOME}" >&2; exit 1; }
mkdir -p "${BUSCO_OUT}" "${PROJECT}/logs"

export CONDA_ENVS_PATH="${PROJECT}/conda_envs"
export CONDA_PKGS_DIRS="${PROJECT}/conda_pkgs"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

busco \
  --in "${PROTEOME}" \
  --out "${SPECIES}_metazoa_odb10" \
  --out_path "${BUSCO_OUT}" \
  --mode proteins \
  --lineage_dataset "${BUSCO_DB}" \
  --cpu "${SLURM_CPUS_PER_TASK}" \
  --offline \
  --force

python3 "${PROJECT}/scripts/06_summarise_primary_busco.py" \
  --busco-dir "${PROJECT}/busco_primary_proteomes" \
  --proteome-dir "${PROJECT}/orthofinder_primary_proteomes" \
  --output "${PROJECT}/primary_proteome_busco_summary.tsv"

grep -E "^${SPECIES}[[:space:]]" "${PROJECT}/primary_proteome_busco_summary.tsv"

