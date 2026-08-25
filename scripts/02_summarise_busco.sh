#!/bin/bash
#SBATCH --job-name=busco_summary
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/busco_summary_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/busco_summary_%j.err

set -euo pipefail

SCRATCH=/path/to/mollusc_project
PROTEOME_DIR="${SCRATCH}/orthofinder_proteomes"
BUSCO_OUT_DIR="${SCRATCH}/busco_proteomes"
SUMMARY_FILE="${SCRATCH}/proteome_busco_summary.tsv"

printf 'species\tproteins\tcomplete_pct\tsingle_copy_pct\tduplicated_pct\tfragmented_pct\tmissing_pct\tbusco_n\n' \
  > "${SUMMARY_FILE}"

for PROTEOME in "${PROTEOME_DIR}"/*.faa; do
  [[ -e "${PROTEOME}" ]] || continue

  SPECIES="$(basename "${PROTEOME}" .faa)"
  N_PROTEINS="$(grep -c '^>' "${PROTEOME}")"
  SUMMARY="$(
    find "${BUSCO_OUT_DIR}/${SPECIES}" -type f -name 'short_summary*.txt' | head -n 1
  )"

  if [[ -z "${SUMMARY}" ]]; then
    printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\n' "${SPECIES}" "${N_PROTEINS}" \
      >> "${SUMMARY_FILE}"
    continue
  fi

  BUSCO_LINE="$(
    grep -E 'C:[0-9.]+%.*S:[0-9.]+%.*D:[0-9.]+%.*F:[0-9.]+%.*M:[0-9.]+%' "${SUMMARY}" \
      | head -n 1 \
      | tr -d '[:space:]'
  )"

  COMPLETE="$(sed -n 's/.*C:\([0-9.]*\)%.*/\1/p' <<< "${BUSCO_LINE}")"
  SINGLE="$(sed -n 's/.*S:\([0-9.]*\)%.*/\1/p' <<< "${BUSCO_LINE}")"
  DUPLICATED="$(sed -n 's/.*D:\([0-9.]*\)%.*/\1/p' <<< "${BUSCO_LINE}")"
  FRAGMENTED="$(sed -n 's/.*F:\([0-9.]*\)%.*/\1/p' <<< "${BUSCO_LINE}")"
  MISSING="$(sed -n 's/.*M:\([0-9.]*\)%.*/\1/p' <<< "${BUSCO_LINE}")"
  BUSCO_N="$(sed -n 's/.*n:\([0-9]*\).*/\1/p' <<< "${BUSCO_LINE}")"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${SPECIES}" "${N_PROTEINS}" "${COMPLETE:-NA}" "${SINGLE:-NA}" \
    "${DUPLICATED:-NA}" "${FRAGMENTED:-NA}" "${MISSING:-NA}" "${BUSCO_N:-NA}" \
    >> "${SUMMARY_FILE}"
done

column -t -s $'\t' "${SUMMARY_FILE}"
echo "BUSCO summary written to ${SUMMARY_FILE}"
