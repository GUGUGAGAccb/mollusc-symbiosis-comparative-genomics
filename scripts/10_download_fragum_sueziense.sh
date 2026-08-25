#!/bin/bash
#SBATCH --job-name=download_fsue
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --account=YOUR_SLURM_ACCOUNT
#SBATCH --output=/path/to/mollusc_project/logs/download_fragum_sueziense_%j.out
#SBATCH --error=/path/to/mollusc_project/logs/download_fragum_sueziense_%j.err

set -euo pipefail

PROJECT=/path/to/mollusc_project
SPECIES=Fragum_sueziense
ACCESSION=GCA_963680895.1
RELEASE=2024_01
SOURCE_DIR="${PROJECT}/source_downloads/${SPECIES}/${ACCESSION}"
GENOME_DIR="${PROJECT}/genome_fastas"
PROTEOME_DIR="${PROJECT}/orthofinder_proteomes"
PRIMARY_DIR="${PROJECT}/orthofinder_primary_proteomes"
MAPPING_DIR="${PROJECT}/primary_isoform_mappings"
TMP_BASE="${PROJECT}/tmp/${SPECIES}_${SLURM_JOB_ID}"

GENOME_GZ="${SPECIES}-${ACCESSION}-unmasked.fa.gz"
PEP_GZ="${SPECIES}-${ACCESSION}-${RELEASE}-pep.fa.gz"
GFF_GZ="${SPECIES}-${ACCESSION}-${RELEASE}-genes.gff3.gz"
GENOME_URL="https://ftp.ensembl.org/pub/rapid-release/species/${SPECIES}/${ACCESSION}/ensembl/genome/${GENOME_GZ}"
GENESET_BASE="https://ftp.ensembl.org/pub/rapid-release/species/${SPECIES}/${ACCESSION}/ensembl/geneset/${RELEASE}"

mkdir -p "${SOURCE_DIR}" "${GENOME_DIR}" "${PROTEOME_DIR}" "${PRIMARY_DIR}" \
  "${MAPPING_DIR}" "${PROJECT}/logs" "${TMP_BASE}/input" "${TMP_BASE}/output" "${TMP_BASE}/mapping"

download() {
  local url="$1"
  local output="$2"
  curl --fail --location --retry 5 --retry-delay 10 --continue-at - \
    --output "${output}" "${url}"
}

download "${GENOME_URL}" "${SOURCE_DIR}/${GENOME_GZ}"
download "${GENESET_BASE}/${PEP_GZ}" "${SOURCE_DIR}/${PEP_GZ}"
download "${GENESET_BASE}/${GFF_GZ}" "${SOURCE_DIR}/${GFF_GZ}"

printf '%s  %s\n' \
  4fb42a5bfb87d5519eb4b1845f790927 "${SOURCE_DIR}/${GENOME_GZ}" \
  78c7c10b0be3722564520ac1d6455b82 "${SOURCE_DIR}/${PEP_GZ}" \
  e21375b5bd0181169a0f77d67634bd24 "${SOURCE_DIR}/${GFF_GZ}" \
  | md5sum --check --strict

gzip -cd "${SOURCE_DIR}/${GENOME_GZ}" > "${TMP_BASE}/${SPECIES}.fna"
gzip -cd "${SOURCE_DIR}/${PEP_GZ}" > "${TMP_BASE}/${SPECIES}.faa"
gzip -cd "${SOURCE_DIR}/${GFF_GZ}" > "${TMP_BASE}/${SPECIES}.gff3"

[[ -s "${TMP_BASE}/${SPECIES}.fna" ]] || { echo "Empty genome FASTA" >&2; exit 1; }
[[ -s "${TMP_BASE}/${SPECIES}.faa" ]] || { echo "Empty protein FASTA" >&2; exit 1; }
[[ -s "${TMP_BASE}/${SPECIES}.gff3" ]] || { echo "Empty GFF3" >&2; exit 1; }

mv "${TMP_BASE}/${SPECIES}.fna" "${GENOME_DIR}/${SPECIES}.fna"
mv "${TMP_BASE}/${SPECIES}.faa" "${PROTEOME_DIR}/${SPECIES}.faa"
mv "${TMP_BASE}/${SPECIES}.gff3" "${SOURCE_DIR}/${SPECIES}.gff3"

ln -s "${PROTEOME_DIR}/${SPECIES}.faa" "${TMP_BASE}/input/${SPECIES}.faa"
"${PROJECT}/conda_envs/orthofinder_env/bin/python3" \
  "${PROJECT}/scripts/04_select_primary_isoforms.py" \
  --input-dir "${TMP_BASE}/input" \
  --output-dir "${TMP_BASE}/output" \
  --mapping-dir "${TMP_BASE}/mapping" \
  --summary "${TMP_BASE}/primary_selection.tsv"

mv "${TMP_BASE}/output/${SPECIES}.faa" "${PRIMARY_DIR}/${SPECIES}.faa"
mv "${TMP_BASE}/mapping/${SPECIES}.removed_isoforms.tsv" \
  "${MAPPING_DIR}/${SPECIES}.removed_isoforms.tsv"
cp "${TMP_BASE}/primary_selection.tsv" "${SOURCE_DIR}/primary_selection.tsv"

cat > "${SOURCE_DIR}/source_metadata.tsv" <<EOF
species\tassembly_accession\tassembly_level\tprovider\tgeneset_release\tgenome_url\tprotein_url\tgff3_url
${SPECIES}\t${ACCESSION}\tchromosome\tEnsembl Rapid Release\t${RELEASE}\t${GENOME_URL}\t${GENESET_BASE}/${PEP_GZ}\t${GENESET_BASE}/${GFF_GZ}
EOF

printf 'Genome sequences: '; grep -c '^>' "${GENOME_DIR}/${SPECIES}.fna"
printf 'Downloaded proteins: '; grep -c '^>' "${PROTEOME_DIR}/${SPECIES}.faa"
printf 'Representative proteins: '; grep -c '^>' "${PRIMARY_DIR}/${SPECIES}.faa"
cat "${TMP_BASE}/primary_selection.tsv"

