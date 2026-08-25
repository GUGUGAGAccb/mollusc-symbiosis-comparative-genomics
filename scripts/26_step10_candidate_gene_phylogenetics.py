#!/usr/bin/env python3
"""Prepare and summarise Step 10 candidate gene-family phylogenetics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median


GROUPS = ("bivalvia", "gastropoda")
STRICT_STATUSES = {"pass", "review"}


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path):
    records = {}
    name = None
    chunks = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper().replace("*", "")
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"Sequence found before FASTA header: {path}")
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks).upper().replace("*", "")
    return records


def number(value, default=0.0):
    try:
        parsed = float(value)
        return parsed if math.isfinite(parsed) else default
    except (TypeError, ValueError):
        return default


def integer(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def safe_label(text):
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_.-")
    return label or "sequence"


def sequence_priority(row, family_median):
    status_rank = 0 if row.get("validation_status") == "pass" else 1
    direct_rank = -integer(row.get("direct_pfam_domain_count"))
    length = integer(row.get("sequence_length_aa"))
    length_distance = abs(length - family_median) if family_median else 0
    return status_rank, direct_rank, length_distance, row.get("query_id", "")


def remove_within_species_duplicates(rows, family_median):
    grouped = defaultdict(list)
    ungrouped = []
    for row in rows:
        duplicate_group = row.get("duplicate_group", "").strip()
        duplicate_context = row.get("duplicate_context", "").strip()
        if duplicate_group and duplicate_context == "within_species_family":
            grouped[(row.get("species", ""), duplicate_group)].append(row)
        else:
            ungrouped.append(row)
    retained = list(ungrouped)
    removed = []
    for members in grouped.values():
        members.sort(key=lambda row: sequence_priority(row, family_median))
        retained.append(members[0])
        removed.extend(members[1:])
    return retained, removed


def cap_family(rows, family_median, per_species_cap, total_cap):
    by_species = defaultdict(list)
    for row in rows:
        by_species[row.get("species", "unknown")].append(row)
    kept = []
    removed = []
    for species in sorted(by_species):
        members = sorted(by_species[species], key=lambda row: sequence_priority(row, family_median))
        kept.extend(members[:per_species_cap])
        removed.extend(members[per_species_cap:])
    kept.sort(key=lambda row: sequence_priority(row, family_median))
    if len(kept) > total_cap:
        removed.extend(kept[total_cap:])
        kept = kept[:total_cap]
    return kept, removed


def family_score(meta, rows, validation_rows):
    evidence = number(meta.get("evidence_score"))
    change = abs(number(meta.get("best_event_change")))
    species = len({row.get("species", "") for row in rows})
    direct = sum(integer(row.get("direct_pfam_domain_count")) > 0 for row in rows)
    review = sum(row.get("validation_status") == "review" for row in rows)
    original = max(1, len(validation_rows))
    direct_fraction = direct / max(1, len(rows))
    review_fraction = review / max(1, len(rows))
    excluded_fraction = (original - len(rows)) / original
    score = (
        evidence
        + min(change, 100.0) * 0.10
        + species * 2.0
        + direct_fraction * 10.0
        - review_fraction * 4.0
        - excluded_fraction * 3.0
    )
    return score, direct_fraction, review_fraction, excluded_fraction


def select_families(eligible, top_n, contraction_quota):
    ordered = sorted(
        eligible,
        key=lambda row: (
            -number(row["selection_score"]),
            -number(row["evidence_score"]),
            -abs(number(row["best_event_change"])),
            row["orthogroup"],
        ),
    )
    contractions = [row for row in ordered if row.get("best_event_direction") == "contraction"]
    selected = contractions[: min(contraction_quota, top_n)]
    selected_ids = {row["orthogroup"] for row in selected}
    for row in ordered:
        if len(selected) >= top_n:
            break
        if row["orthogroup"] not in selected_ids:
            selected.append(row)
            selected_ids.add(row["orthogroup"])
    selected.sort(
        key=lambda row: (
            -number(row["selection_score"]),
            -number(row["evidence_score"]),
            row["orthogroup"],
        )
    )
    for rank, row in enumerate(selected, 1):
        row["rank"] = rank
    return selected


def prepare_group(args, group, run_dir):
    validation_dir = args.project_root / "results_or_reports" / "candidate_sequence_validation" / group / "pfam_hmmer_v1" / "latest"
    annotation_dir = args.project_root / "results_or_reports" / "functional_annotation" / group / "eggnog_mapper_v2_1_12" / "latest"
    validation_path = validation_dir / "sequence_validation.tsv"
    family_path = annotation_dir / "high_confidence_family_functional_annotations.tsv"
    fasta_path = annotation_dir / "annotation_input_unique.faa"
    for path in (validation_dir / "COMPLETE", validation_path, family_path, fasta_path):
        if not path.exists() or (path.is_file() and path.stat().st_size == 0):
            raise FileNotFoundError(f"Required Step 9/7 input missing: {path}")

    validation = read_tsv(validation_path)
    family_meta = {row["orthogroup"]: row for row in read_tsv(family_path)}
    sequences = read_fasta(fasta_path)
    by_family = defaultdict(list)
    for row in validation:
        by_family[row["orthogroup"]].append(row)

    eligible = []
    retained_by_family = {}
    eligibility_rows = []
    for orthogroup, meta in sorted(family_meta.items()):
        all_rows = by_family.get(orthogroup, [])
        strict = [row for row in all_rows if row.get("validation_status") in STRICT_STATUSES]
        lengths = [integer(row.get("sequence_length_aa")) for row in strict if integer(row.get("sequence_length_aa")) > 0]
        family_median = median(lengths) if lengths else 0
        strict, duplicate_removed = remove_within_species_duplicates(strict, family_median)
        strict, cap_removed = cap_family(strict, family_median, args.per_species_cap, args.total_cap)
        strict = [row for row in strict if row.get("query_id") in sequences and len(sequences[row["query_id"]]) >= args.min_length]
        species_count = len({row.get("species", "") for row in strict})
        unique_sequence_count = len({hashlib.md5(sequences[row["query_id"]].encode()).hexdigest() for row in strict})
        reasons = []
        if len(strict) < args.min_sequences:
            reasons.append("too_few_strict_sequences")
        if species_count < args.min_species:
            reasons.append("too_few_species")
        if unique_sequence_count < args.min_unique_sequences:
            reasons.append("too_few_unique_sequences")
        eligible_flag = not reasons
        score, direct_fraction, review_fraction, excluded_fraction = family_score(meta, strict, all_rows)
        row = {
            "group": group,
            "orthogroup": orthogroup,
            "eligible": "yes" if eligible_flag else "no",
            "ineligibility_reason": ";".join(reasons),
            "strict_sequences": len(strict),
            "strict_species": species_count,
            "unique_sequences": unique_sequence_count,
            "original_sequences": len(all_rows),
            "within_species_duplicates_removed": len(duplicate_removed),
            "cap_removed": len(cap_removed),
            "direct_pfam_fraction": f"{direct_fraction:.6f}",
            "review_fraction": f"{review_fraction:.6f}",
            "excluded_fraction": f"{excluded_fraction:.6f}",
            "selection_score": f"{score:.6f}",
            "evidence_score": meta.get("evidence_score", ""),
            "candidate_tier": meta.get("candidate_tier", ""),
            "distribution_category": meta.get("distribution_category", ""),
            "best_event_node": meta.get("best_event_node", ""),
            "best_event_direction": meta.get("best_event_direction", ""),
            "best_event_change": meta.get("best_event_change", ""),
            "family_q_bh": meta.get("family_q_bh", ""),
            "consensus_preferred_name": meta.get("consensus_preferred_name", ""),
            "consensus_description": meta.get("consensus_description", ""),
            "pfam_domains": meta.get("pfam_domains", ""),
        }
        eligibility_rows.append(row)
        if eligible_flag:
            eligible.append(row)
            retained_by_family[orthogroup] = strict

    selected = select_families(eligible, args.top_n, args.contraction_quota)
    input_dir = run_dir / "inputs" / group
    input_dir.mkdir(parents=True, exist_ok=True)
    selected_rows = []
    for row in selected:
        orthogroup = row["orthogroup"]
        members = retained_by_family[orthogroup]
        fasta_out = input_dir / f"{orthogroup}.faa"
        map_out = input_dir / f"{orthogroup}.sequence_map.tsv"
        label_counts = Counter()
        map_rows = []
        with fasta_out.open("w") as handle:
            for member in sorted(members, key=lambda item: (item.get("species", ""), item.get("gene_id", ""), item.get("query_id", ""))):
                base = safe_label(f"{member.get('species', 'unknown')}__{member.get('gene_id', member['query_id'])}")
                label_counts[base] += 1
                label = base if label_counts[base] == 1 else f"{base}__{label_counts[base]}"
                sequence = sequences[member["query_id"]]
                handle.write(f">{label}\n")
                for start in range(0, len(sequence), 60):
                    handle.write(sequence[start:start + 60] + "\n")
                map_rows.append({
                    "tree_label": label,
                    "query_id": member["query_id"],
                    "species": member.get("species", ""),
                    "gene_id": member.get("gene_id", ""),
                    "orthogroup": orthogroup,
                    "validation_status": member.get("validation_status", ""),
                    "sequence_length_aa": len(sequence),
                    "direct_pfam_accessions": member.get("direct_pfam_accessions", ""),
                })
        write_tsv(map_out, ["tree_label", "query_id", "species", "gene_id", "orthogroup", "validation_status", "sequence_length_aa", "direct_pfam_accessions"], map_rows)
        selected_rows.append({
            "rank": row["rank"],
            "group": group,
            "orthogroup": orthogroup,
            "input_fasta": str(fasta_out),
            "sequence_map": str(map_out),
            "n_sequences": len(map_rows),
            "n_species": len({item["species"] for item in map_rows}),
            **{key: value for key, value in row.items() if key not in {"rank", "group", "orthogroup"}},
        })

    eligibility_fields = list(eligibility_rows[0]) if eligibility_rows else ["group", "orthogroup", "eligible"]
    selected_fields = list(selected_rows[0]) if selected_rows else ["rank", "group", "orthogroup", "input_fasta", "sequence_map", "n_sequences", "n_species"]
    write_tsv(run_dir / f"eligible_{group}.tsv", eligibility_fields, eligibility_rows)
    write_tsv(run_dir / f"selected_{group}.tsv", selected_fields, selected_rows)
    return {
        "group": group,
        "input_candidate_families": len(family_meta),
        "eligible_families": len(eligible),
        "selected_families": len(selected),
        "selected_expansions": sum(row.get("best_event_direction") == "expansion" for row in selected),
        "selected_contractions": sum(row.get("best_event_direction") == "contraction" for row in selected),
        "selected_sequences": sum(integer(row.get("n_sequences")) for row in selected_rows),
        "selected_species_min": min((integer(row.get("n_species")) for row in selected_rows), default=0),
        "selected_species_max": max((integer(row.get("n_species")) for row in selected_rows), default=0),
    }


def cmd_prepare(args):
    args.run_dir.mkdir(parents=True, exist_ok=True)
    summaries = [prepare_group(args, group, args.run_dir) for group in GROUPS]
    write_tsv(args.run_dir / "selection_summary.tsv", list(summaries[0]), summaries)
    methods = f"""# Step 10 candidate gene-family phylogenetics

Candidate families were drawn independently from the Bivalvia and Gastropoda high-confidence family tables.
Only Step 9 sequences classified as `pass` or ordinary `review` were eligible. `review_high` and
`exclude_recommended` sequences were excluded. Exact duplicates within the same species and family were
collapsed, retaining the better validated representative. Families required at least {args.min_sequences}
sequences, {args.min_species} species and {args.min_unique_sequences} unique amino-acid sequences.

To bound computational cost and reduce domination by extreme expansions, at most {args.per_species_cap}
sequences per species and {args.total_cap} sequences per family were retained using validation status,
direct-Pfam support and proximity to family median length. The top {args.top_n} eligible families per class
were ranked by family-evolution evidence, event magnitude, taxon representation, direct-Pfam coverage and
sequence-validation penalties. Up to {args.contraction_quota} places were reserved for contraction families.

Downstream tree tasks use MAFFT, trimAl and IQ-TREE restricted ModelFinder with 1,000 ultrafast bootstrap
and 1,000 SH-aLRT replicates. If support estimation fails for a very difficult family, an ML-only fallback
is attempted and explicitly recorded. These gene trees are exploratory and require Step 11 reconciliation
against the species tree before duplication/loss interpretation.
"""
    (args.run_dir / "METHODS_AND_LIMITATIONS.txt").write_text(methods)
    (args.run_dir / "PREPARED").write_text("Step 10 inputs prepared\n")


def parse_iqtree(path: Path):
    result = {"best_model": "", "log_likelihood": "", "alignment_sequences": "", "alignment_sites": "", "cpu_time": "", "wall_time": ""}
    if not path.exists():
        return result
    text = path.read_text(errors="replace")
    patterns = {
        "best_model": r"Best-fit model according to BIC:\s*(\S+)",
        "log_likelihood": r"Log-likelihood of the tree:\s*([-0-9.eE+]+)",
        "alignment_sequences": r"Alignment has\s+(\d+)\s+sequences",
        "alignment_sites": r"Alignment has\s+\d+\s+sequences with\s+(\d+)\s+columns",
        "cpu_time": r"Total CPU time used:\s*(.+)",
        "wall_time": r"Total wall-clock time used:\s*(.+)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if match:
            result[key] = match.group(1).strip()
    return result


def cmd_summarise(args):
    run_dir = args.run_dir.resolve()
    tree_rows = []
    failed_rows = []
    for group in GROUPS:
        manifest = read_tsv(run_dir / f"selected_{group}.tsv")
        for selected in manifest:
            orthogroup = selected["orthogroup"]
            tree_dir = run_dir / "trees" / group / orthogroup
            prefix_file = tree_dir / "TREE_PREFIX"
            status = "missing"
            support_mode = ""
            prefix = ""
            if (tree_dir / "COMPLETE").exists() and prefix_file.exists():
                prefix = prefix_file.read_text().strip()
                support_mode = (tree_dir / "SUPPORT_MODE").read_text().strip() if (tree_dir / "SUPPORT_MODE").exists() else "unknown"
                status = "complete"
            elif (tree_dir / "FAILED").exists():
                status = "failed"
            iqtree_path = Path(prefix + ".iqtree") if prefix else Path("/__missing__")
            treefile = Path(prefix + ".treefile") if prefix else Path("/__missing__")
            contree = Path(prefix + ".contree") if prefix else Path("/__missing__")
            parsed = parse_iqtree(iqtree_path)
            row = {
                "group": group,
                "rank": selected.get("rank", ""),
                "orthogroup": orthogroup,
                "status": status,
                "support_mode": support_mode,
                "input_sequences": selected.get("n_sequences", ""),
                "input_species": selected.get("n_species", ""),
                "best_event_direction": selected.get("best_event_direction", ""),
                "best_event_change": selected.get("best_event_change", ""),
                "selection_score": selected.get("selection_score", ""),
                "best_model": parsed["best_model"],
                "log_likelihood": parsed["log_likelihood"],
                "alignment_sequences": parsed["alignment_sequences"],
                "alignment_sites": parsed["alignment_sites"],
                "cpu_time": parsed["cpu_time"],
                "wall_time": parsed["wall_time"],
                "treefile": str(treefile) if treefile.exists() else "",
                "contree": str(contree) if contree.exists() else "",
                "iqtree_report": str(iqtree_path) if iqtree_path.exists() else "",
                "sequence_map": selected.get("sequence_map", ""),
                "failure_message": (tree_dir / "FAILED").read_text(errors="replace").strip() if (tree_dir / "FAILED").exists() else "",
            }
            tree_rows.append(row)
            if status != "complete" or not treefile.exists():
                failed_rows.append(row)
    fields = list(tree_rows[0]) if tree_rows else ["group", "orthogroup", "status"]
    write_tsv(run_dir / "gene_tree_summary.tsv", fields, tree_rows)
    write_tsv(run_dir / "failed_or_missing_trees.tsv", fields, failed_rows)
    summary_rows = []
    for group in GROUPS:
        rows = [row for row in tree_rows if row["group"] == group]
        summary_rows.append({
            "group": group,
            "selected_families": len(rows),
            "completed_trees": sum(row["status"] == "complete" for row in rows),
            "supported_trees": sum(row["support_mode"] == "ufboot1000_alrt1000" for row in rows),
            "ml_only_fallback_trees": sum(row["support_mode"] == "ml_only_fallback" for row in rows),
            "failed_or_missing_trees": sum(row["status"] != "complete" for row in rows),
        })
    write_tsv(run_dir / "step10_summary.tsv", list(summary_rows[0]), summary_rows)
    manifest_rows = []
    for path in sorted(run_dir.rglob("*")):
        if path.is_file() and path.name not in {"result_manifest.tsv"}:
            manifest_rows.append({"relative_path": str(path.relative_to(run_dir)), "bytes": path.stat().st_size})
    write_tsv(run_dir / "result_manifest.tsv", ["relative_path", "bytes"], manifest_rows)
    (run_dir / "COMPLETE").write_text(f"Step 10 summary complete; failed_or_missing={len(failed_rows)}\n")


def build_parser():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare")
    prepare.add_argument("--project-root", type=Path, required=True)
    prepare.add_argument("--run-dir", type=Path, required=True)
    prepare.add_argument("--top-n", type=int, default=30)
    prepare.add_argument("--contraction-quota", type=int, default=6)
    prepare.add_argument("--min-sequences", type=int, default=4)
    prepare.add_argument("--min-species", type=int, default=3)
    prepare.add_argument("--min-unique-sequences", type=int, default=3)
    prepare.add_argument("--min-length", type=int, default=30)
    prepare.add_argument("--per-species-cap", type=int, default=20)
    prepare.add_argument("--total-cap", type=int, default=180)
    prepare.set_defaults(func=cmd_prepare)
    summarise = sub.add_parser("summarise")
    summarise.add_argument("--run-dir", type=Path, required=True)
    summarise.set_defaults(func=cmd_summarise)
    return parser


if __name__ == "__main__":
    parsed = build_parser().parse_args()
    parsed.func(parsed)
