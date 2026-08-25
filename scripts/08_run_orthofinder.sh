#!/bin/bash
#SBATCH --job-name=orthofinder_core
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=3-00:00:00
#SBATCH --mem=80G
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/orthofinder_core_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/orthofinder_core_%j.err

set -euo pipefail

SCRATCH=/path/to/mollusc_project
PRIMARY_DIR="${SCRATCH}/orthofinder_primary_proteomes"
CORE_DIR="${SCRATCH}/orthofinder_core_complete80"
RESULTS_PARENT="${SCRATCH}/orthofinder_results"
RESULTS_DIR="${RESULTS_PARENT}/core_complete80_v3_1_5"
SUMMARY="${SCRATCH}/primary_proteome_busco_summary.tsv"
DECISIONS="${SCRATCH}/orthofinder_species_decisions.tsv"
CONDA_BASE=/path/to/miniconda3
CONDA_ENV="${SCRATCH}/conda_envs/orthofinder_env"

export CONDA_ENVS_PATH="${SCRATCH}/conda_envs"
export CONDA_PKGS_DIRS="${SCRATCH}/conda_pkgs"
export TMPDIR="${SCRATCH}/tmp/orthofinder_${SLURM_JOB_ID}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

mkdir -p "${CORE_DIR}" "${RESULTS_PARENT}" "${TMPDIR}"

python3 "${SCRATCH}/scripts/07_prepare_core_set.py" \
  --summary "${SUMMARY}" \
  --proteome-dir "${PRIMARY_DIR}" \
  --core-dir "${CORE_DIR}" \
  --decision-table "${DECISIONS}" \
  --min-complete 80

if [[ -e "${RESULTS_DIR}" ]]; then
  echo "ERROR: results path already exists: ${RESULTS_DIR}" >&2
  exit 1
fi

column -t -s $'\t' "${DECISIONS}"
orthofinder \
  -f "${CORE_DIR}" \
  -t "${SLURM_CPUS_PER_TASK}" \
  -a "${SLURM_CPUS_PER_TASK}" \
  -M msa \
  -S diamond \
  -A famsa \
  -T fasttree \
  -o "${RESULTS_DIR}"
