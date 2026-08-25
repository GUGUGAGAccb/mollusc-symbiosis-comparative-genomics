#!/bin/bash
#SBATCH --job-name=primary_busco
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --array=1-14%4
#SBATCH --output=/path/to/mollusc_project/logs/primary_busco_%A_%a.out
#SBATCH --error=/path/to/mollusc_project/logs/primary_busco_%A_%a.err

set -euo pipefail

SCRATCH=/path/to/mollusc_project
PROTEOME_DIR="${SCRATCH}/orthofinder_primary_proteomes"
BUSCO_DB="${SCRATCH}/busco_downloads/lineages/metazoa_odb10"
BUSCO_OUT="${SCRATCH}/busco_primary_proteomes"
CONDA_BASE=/path/to/miniconda3
CONDA_ENV="${SCRATCH}/conda_envs/busco_env"

export CONDA_ENVS_PATH="${SCRATCH}/conda_envs"
export CONDA_PKGS_DIRS="${SCRATCH}/conda_pkgs"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

mapfile -t PROTEOMES < <(find "${PROTEOME_DIR}" -maxdepth 1 -type f -name '*.faa' -size +0c | sort)
TASK_ID="${SLURM_ARRAY_TASK_ID}"
(( TASK_ID >= 1 && TASK_ID <= ${#PROTEOMES[@]} )) || { echo "Invalid array task ${TASK_ID}" >&2; exit 1; }

PROTEOME="${PROTEOMES[$((TASK_ID - 1))]}"
SPECIES="$(basename "${PROTEOME}" .faa)"
SPECIES_OUT="${BUSCO_OUT}/${SPECIES}"
mkdir -p "${SPECIES_OUT}"

busco \
  --in "${PROTEOME}" \
  --out "${SPECIES}_metazoa_odb10" \
  --out_path "${SPECIES_OUT}" \
  --mode proteins \
  --lineage_dataset "${BUSCO_DB}" \
  --cpu "${SLURM_CPUS_PER_TASK}" \
  --offline \
  --force
