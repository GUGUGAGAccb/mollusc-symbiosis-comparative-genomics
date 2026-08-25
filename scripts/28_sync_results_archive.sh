#!/usr/bin/env bash
set -euo pipefail

ROOT=${MOLLUSC_PROJECT_ROOT:-/path/to/mollusc_project}
RESULT_ROOT=${RESULT_ROOT:-$ROOT/result}
MODE=${1:-all}
mkdir -p "$RESULT_ROOT"

exec 9>"$RESULT_ROOT/.archive.lock"
flock -w 3600 9

SOURCE_TABLE="$RESULT_ROOT/archive_sources.tsv"
if [[ "$MODE" == all ]]; then
  printf 'step\tlabel\tsource\tresolved_source\tdestination\tstatus\n' > "$SOURCE_TABLE"
elif [[ ! -s "$SOURCE_TABLE" ]]; then
  printf 'step\tlabel\tsource\tresolved_source\tdestination\tstatus\n' > "$SOURCE_TABLE"
fi

record_source() {
  local step=$1 label=$2 source=$3 resolved=$4 destination=$5 status=$6
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$step" "$label" "$source" "$resolved" "$destination" "$status" >> "$SOURCE_TABLE"
}

archive_dir() {
  local source=$1 step=$2 label=$3
  shift 3
  local resolved destination
  if [[ ! -e "$source" ]]; then
    record_source "$step" "$label" "$source" "" "" "missing"
    printf 'WARNING: source missing: %s\n' "$source" >&2
    return 0
  fi
  resolved=$(readlink -f "$source")
  destination="$RESULT_ROOT/$step/$label"
  mkdir -p "$destination"
  rsync -aL --partial "$@" "$resolved/" "$destination/"
  ln -sfn "$label" "$RESULT_ROOT/$step/latest"
  record_source "$step" "$label" "$source" "$resolved" "$destination" "copied"
}

archive_file() {
  local source=$1 step=$2 label=${3:-current}
  local resolved destination
  if [[ ! -f "$source" ]]; then
    record_source "$step" "$label" "$source" "" "" "missing"
    return 0
  fi
  resolved=$(readlink -f "$source")
  destination="$RESULT_ROOT/$step/$label"
  mkdir -p "$destination"
  rsync -a "$resolved" "$destination/"
  record_source "$step" "$label" "$source" "$resolved" "$destination" "copied"
}

archive_tar_gz() {
  local source=$1 step=$2 label=$3 archive_name=$4
  local resolved destination temporary
  if [[ ! -d "$source" ]]; then
    record_source "$step" "$label" "$source" "" "" "missing"
    return 0
  fi
  resolved=$(readlink -f "$source")
  destination="$RESULT_ROOT/$step/$label"
  mkdir -p "$destination"
  temporary="$destination/${archive_name}.partial"
  tar -C "$(dirname "$resolved")" -czf "$temporary" "$(basename "$resolved")"
  mv -f "$temporary" "$destination/$archive_name"
  record_source "$step" "$label" "$source" "$resolved" "$destination/$archive_name" "copied_compressed"
}

write_metadata() {
  cat > "$RESULT_ROOT/ARCHIVE_POLICY.md" <<'EOF'
# Mollusc project central results archive

This directory contains copied final results organised by pipeline step. Original analysis directories remain
unchanged and authoritative. The archive excludes raw downloaded genomes/proteomes, databases, Conda environments,
failed retry directories, Slurm logs and OrthoFinder `WorkingDirectory` caches. Completed final tables, trees,
candidate sequences, annotations, enrichment outputs and method/limitation files are copied. To remain within the
BluePebble quota, Step 02 keeps final OrthoFinder tables, statistics, species trees and gene-duplication outputs;
its large `Gene_Duplication_Events` and `Phylogenetic_Hierarchical_Orthogroups` result directories are retained as
lossless `.tar.gz` copies. The reproducible per-orthogroup sequence, resolved-gene-tree and pairwise-orthologue
collections are not duplicated and remain unchanged in the authoritative OrthoFinder source directories.

The archive is append-safe and idempotent: repeated synchronization does not delete prior run directories. Each
step or subgroup has a `latest` symlink pointing to the most recently synchronized run. Future steps should call:

`28_sync_results_archive.sh one stepNN_descriptive_name /absolute/path/to/completed/run run_JOBID`

Only completed outputs should be passed to the `one` mode. `archive_sources.tsv` records provenance and
`archive_inventory.tsv` records file size and modification time for verification.
EOF
  date -Is > "$RESULT_ROOT/LAST_SYNC_TIMESTAMP.txt"
  python - "$RESULT_ROOT" <<'PY'
import os
import sys
from pathlib import Path

root = Path(sys.argv[1])
inventory = root / "archive_inventory.tsv"
summary = root / "archive_step_summary.tsv"
excluded = {inventory.name, summary.name, ".archive.lock"}
rows = []
step_counts = {}
for path in sorted(root.rglob("*")):
    if not path.is_file() or path.name in excluded:
        continue
    stat = path.stat()
    rel = path.relative_to(root)
    rows.append((str(rel), stat.st_size, stat.st_mtime_ns))
    step = rel.parts[0] if rel.parts else "root_metadata"
    count, size = step_counts.get(step, (0, 0))
    step_counts[step] = (count + 1, size + stat.st_size)
# Some historical analysis files can contain undecodable byte sequences in
# their names. Preserve a reversible escaped representation in the TSV rather
# than allowing one legacy filename to abort the entire archive.
with inventory.open("w", encoding="utf-8", errors="backslashreplace") as handle:
    handle.write("relative_path\tbytes\tmtime_ns\n")
    for rel, size, mtime in rows:
        handle.write(f"{rel}\t{size}\t{mtime}\n")
with summary.open("w", encoding="utf-8", errors="backslashreplace") as handle:
    handle.write("step\tfiles\tbytes\tmegabytes\n")
    for step, (count, size) in sorted(step_counts.items()):
        handle.write(f"{step}\t{count}\t{size}\t{size / 1048576:.3f}\n")
PY
  printf 'Central result archive synchronized successfully.\n' > "$RESULT_ROOT/COPY_COMPLETE"
}

if [[ "$MODE" == one ]]; then
  if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: $0 one STEP_NAME SOURCE_DIR [RUN_LABEL]" >&2
    exit 2
  fi
  STEP_NAME=$2
  SOURCE_DIR=$3
  RUN_LABEL=${4:-$(basename "$(readlink -f "$SOURCE_DIR")")}
  archive_dir "$SOURCE_DIR" "$STEP_NAME" "$RUN_LABEL"
  write_metadata
  exit 0
fi

if [[ "$MODE" != all ]]; then
  echo "Usage: $0 [all] | one STEP_NAME SOURCE_DIR [RUN_LABEL]" >&2
  exit 2
fi

# A successful marker must describe the current full synchronization, not an
# earlier preflight copy. Preserve the previous marker for auditability while
# the new synchronization is running.
if [[ -f "$RESULT_ROOT/COPY_COMPLETE" ]]; then
  mv -f "$RESULT_ROOT/COPY_COMPLETE" "$RESULT_ROOT/COPY_COMPLETE.previous"
fi
printf 'RUNNING\t%s\n' "$(date -Is)" > "$RESULT_ROOT/ARCHIVE_STATUS.tsv"
archive_finished=0
archive_exit_status() {
  local rc=$?
  if (( archive_finished == 0 )); then
    printf 'FAILED\t%s\texit_code=%s\n' "$(date -Is)" "$rc" > "$RESULT_ROOT/ARCHIVE_STATUS.tsv" 2>/dev/null || true
  fi
}
trap archive_exit_status EXIT

# Step 01: proteome quality, representative-proteome selection and inclusion decisions.
for file in \
  genome_proteome_manifest.tsv \
  orthofinder_species_decisions.tsv \
  primary_proteome_busco_summary.tsv \
  primary_proteome_selection.tsv \
  proteome_busco_summary.tsv; do
  archive_file "$ROOT/$file" "step01_proteome_quality" "current"
done

# Step 02: quota-safe final OrthoFinder deliverables. Large reproducible
# sequence/tree collections remain in the authoritative source directories.
archive_dir "$ROOT/orthofinder_results/core_complete80_v3_1_5/Results_Aug04" \
  "step02_orthogroup_inference/combined_14_species" "Results_Aug04" \
  --exclude='/WorkingDirectory/' --exclude='/Orthologues/' --exclude='/Orthogroup_Sequences/' \
  --exclude='/Resolved_Gene_Trees/' --exclude='/Single_Copy_Orthologue_Sequences/' \
  --exclude='/MultipleSequenceAlignments/' --exclude='/Gene_Duplication_Events/' \
  --exclude='/Phylogenetic_Hierarchical_Orthogroups/'
archive_dir "$ROOT/orthofinder_results/by_class/bivalvia_v3_1_5/Results_Aug04" \
  "step02_orthogroup_inference/bivalvia" "Results_Aug04" \
  --exclude='/WorkingDirectory/' --exclude='/Orthologues/' --exclude='/Orthogroup_Sequences/' \
  --exclude='/Resolved_Gene_Trees/' --exclude='/Single_Copy_Orthologue_Sequences/' \
  --exclude='/MultipleSequenceAlignments/' --exclude='/Gene_Duplication_Events/' \
  --exclude='/Phylogenetic_Hierarchical_Orthogroups/'
archive_dir "$ROOT/orthofinder_results/by_class/gastropoda_v3_1_5/Results_Aug04" \
  "step02_orthogroup_inference/gastropoda" "Results_Aug04" \
  --exclude='/WorkingDirectory/' --exclude='/Orthologues/' --exclude='/Orthogroup_Sequences/' \
  --exclude='/Resolved_Gene_Trees/' --exclude='/Single_Copy_Orthologue_Sequences/' \
  --exclude='/MultipleSequenceAlignments/' --exclude='/Gene_Duplication_Events/' \
  --exclude='/Phylogenetic_Hierarchical_Orthogroups/'

for group in combined_14_species bivalvia gastropoda; do
  case "$group" in
    combined_14_species) of_root="$ROOT/orthofinder_results/core_complete80_v3_1_5/Results_Aug04" ;;
    bivalvia) of_root="$ROOT/orthofinder_results/by_class/bivalvia_v3_1_5/Results_Aug04" ;;
    gastropoda) of_root="$ROOT/orthofinder_results/by_class/gastropoda_v3_1_5/Results_Aug04" ;;
  esac
  archive_tar_gz "$of_root/Gene_Duplication_Events" \
    "step02_orthogroup_inference/$group" "Results_Aug04" "Gene_Duplication_Events.tar.gz"
  archive_tar_gz "$of_root/Phylogenetic_Hierarchical_Orthogroups" \
    "step02_orthogroup_inference/$group" "Results_Aug04" "Phylogenetic_Hierarchical_Orthogroups.tar.gz"
done

# Step 03: explicit species-tree copies organised separately from OrthoFinder results.
archive_dir "$ROOT/orthofinder_results/core_complete80_v3_1_5/Results_Aug04/Species_Tree" \
  "step03_species_tree/combined_14_species/orthofinder_rooted" "Results_Aug04"
archive_dir "$ROOT/orthofinder_results/core_complete80_v3_1_5/Results_Aug04/Species_Tree_IQTree_restricted_mfp_v1" \
  "step03_species_tree/combined_14_species/iqtree_restricted" "Results_Aug04"
archive_dir "$ROOT/orthofinder_results/by_class/bivalvia_v3_1_5/Results_Aug04/Species_Tree" \
  "step03_species_tree/bivalvia/orthofinder_rooted" "Results_Aug04"
archive_dir "$ROOT/orthofinder_results/by_class/bivalvia_v3_1_5/Results_Aug04/Species_Tree_IQTree" \
  "step03_species_tree/bivalvia/iqtree" "Results_Aug04"
archive_dir "$ROOT/orthofinder_results/by_class/gastropoda_v3_1_5/Results_Aug04/Species_Tree" \
  "step03_species_tree/gastropoda/orthofinder_rooted" "Results_Aug04"
archive_dir "$ROOT/orthofinder_results/by_class/gastropoda_v3_1_5/Results_Aug04/Species_Tree_IQTree_restricted_mfp_v1" \
  "step03_species_tree/gastropoda/iqtree_restricted" "Results_Aug04"

# Steps 04-11: copy the resolved successful/latest run directories in full.
archive_dir "$ROOT/results_or_reports/orthogroup_distribution/latest" \
  "step04_orthogroup_distribution" "$(basename "$(readlink -f "$ROOT/results_or_reports/orthogroup_distribution/latest")")"
archive_dir "$ROOT/results_or_reports/gene_family_evolution/bivalvia/relative_time_screen_v1/latest" \
  "step05_gene_family_evolution/bivalvia" "$(basename "$(readlink -f "$ROOT/results_or_reports/gene_family_evolution/bivalvia/relative_time_screen_v1/latest")")"
archive_dir "$ROOT/results_or_reports/gene_family_evolution/gastropoda/relative_time_screen_v1/latest" \
  "step05_gene_family_evolution/gastropoda" "$(basename "$(readlink -f "$ROOT/results_or_reports/gene_family_evolution/gastropoda/relative_time_screen_v1/latest")")"
archive_dir "$ROOT/results_or_reports/candidate_gene_families/base_model_v1/latest" \
  "step06_candidate_gene_families" "$(basename "$(readlink -f "$ROOT/results_or_reports/candidate_gene_families/base_model_v1/latest")")"
archive_dir "$ROOT/results_or_reports/functional_annotation/bivalvia/eggnog_mapper_v2_1_12/latest" \
  "step07_functional_annotation/bivalvia" "$(basename "$(readlink -f "$ROOT/results_or_reports/functional_annotation/bivalvia/eggnog_mapper_v2_1_12/latest")")"
archive_dir "$ROOT/results_or_reports/functional_annotation/gastropoda/eggnog_mapper_v2_1_12/latest" \
  "step07_functional_annotation/gastropoda" "$(basename "$(readlink -f "$ROOT/results_or_reports/functional_annotation/gastropoda/eggnog_mapper_v2_1_12/latest")")"
archive_dir "$ROOT/results_or_reports/functional_annotation/combined_summary/latest" \
  "step07_functional_annotation/combined" "$(basename "$(readlink -f "$ROOT/results_or_reports/functional_annotation/combined_summary/latest")")"
archive_dir "$ROOT/results_or_reports/functional_enrichment/orthogroup_level_v1/latest" \
  "step08_functional_enrichment" "$(basename "$(readlink -f "$ROOT/results_or_reports/functional_enrichment/orthogroup_level_v1/latest")")"
archive_dir "$ROOT/results_or_reports/candidate_sequence_validation/bivalvia/pfam_hmmer_v1/latest" \
  "step09_candidate_sequence_validation/bivalvia" "$(basename "$(readlink -f "$ROOT/results_or_reports/candidate_sequence_validation/bivalvia/pfam_hmmer_v1/latest")")"
archive_dir "$ROOT/results_or_reports/candidate_sequence_validation/gastropoda/pfam_hmmer_v1/latest" \
  "step09_candidate_sequence_validation/gastropoda" "$(basename "$(readlink -f "$ROOT/results_or_reports/candidate_sequence_validation/gastropoda/pfam_hmmer_v1/latest")")"
archive_dir "$ROOT/results_or_reports/candidate_sequence_validation/combined_summary/latest" \
  "step09_candidate_sequence_validation/combined" "$(basename "$(readlink -f "$ROOT/results_or_reports/candidate_sequence_validation/combined_summary/latest")")"
archive_dir "$ROOT/results_or_reports/candidate_family_phylogenetics/validated_top30_v1/latest" \
  "step10_candidate_family_phylogenetics" "$(basename "$(readlink -f "$ROOT/results_or_reports/candidate_family_phylogenetics/validated_top30_v1/latest")")"
archive_dir "$ROOT/results_or_reports/gene_tree_species_tree_reconciliation/bivalvia/dl_lca_v1/latest" \
  "step11_gene_tree_species_tree_reconciliation/bivalvia" "$(basename "$(readlink -f "$ROOT/results_or_reports/gene_tree_species_tree_reconciliation/bivalvia/dl_lca_v1/latest")")"
archive_dir "$ROOT/results_or_reports/gene_tree_species_tree_reconciliation/gastropoda/dl_lca_v1/latest" \
  "step11_gene_tree_species_tree_reconciliation/gastropoda" "$(basename "$(readlink -f "$ROOT/results_or_reports/gene_tree_species_tree_reconciliation/gastropoda/dl_lca_v1/latest")")"
archive_dir "$ROOT/results_or_reports/gene_tree_species_tree_reconciliation/combined_summary/dl_lca_v1/latest" \
  "step11_gene_tree_species_tree_reconciliation/combined" "$(basename "$(readlink -f "$ROOT/results_or_reports/gene_tree_species_tree_reconciliation/combined_summary/dl_lca_v1/latest")")"

# Copy the pipeline scripts required to reproduce or extend the archived results.
mkdir -p "$RESULT_ROOT/reproducibility/scripts"
rsync -a "$ROOT/scripts/" "$RESULT_ROOT/reproducibility/scripts/"
record_source "reproducibility" "scripts" "$ROOT/scripts" "$ROOT/scripts" "$RESULT_ROOT/reproducibility/scripts" "copied"

write_metadata
archive_finished=1
printf 'COMPLETED\t%s\n' "$(date -Is)" > "$RESULT_ROOT/ARCHIVE_STATUS.tsv"
printf 'All completed pipeline results synchronized to %s\n' "$RESULT_ROOT"
