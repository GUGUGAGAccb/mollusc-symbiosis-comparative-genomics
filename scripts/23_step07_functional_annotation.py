#!/usr/bin/env python3
"""Prepare eggNOG-mapper input and summarise Step 7 functional annotations."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path


MISSING = {"", "-", "NA", "N/A", "None", "none", "null"}


def read_fasta(path: Path):
    name = None
    chunks = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def parse_query_id(query_id: str, expected_group: str):
    parts = query_id.split("|")
    if len(parts) != 3:
        raise ValueError(f"Expected species|gene|orthogroup query ID, found: {query_id}")
    species, gene_id, orthogroup = parts
    if not species or not gene_id or not orthogroup.startswith("OG"):
        raise ValueError(f"Malformed query ID for {expected_group}: {query_id}")
    return species, gene_id, orthogroup


def write_tsv(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normal(value):
    value = (value or "").strip()
    return "" if value in MISSING else value


def split_terms(value):
    value = normal(value)
    if not value:
        return []
    return [term.strip() for term in value.split(",") if term.strip() and term.strip() not in MISSING]


def cmd_prepare(args):
    group = args.group
    records = OrderedDict()

    def add(path: Path, source: str):
        for query_id, sequence in read_fasta(path):
            species, gene_id, orthogroup = parse_query_id(query_id, group)
            if not sequence:
                raise ValueError(f"Empty protein sequence: {query_id}")
            if query_id not in records:
                records[query_id] = {
                    "sequence": sequence,
                    "species": species,
                    "gene_id": gene_id,
                    "orthogroup": orthogroup,
                    "sources": set(),
                }
            elif records[query_id]["sequence"] != sequence:
                raise ValueError(f"Conflicting sequences for query ID: {query_id}")
            records[query_id]["sources"].add(source)

    add(args.high_confidence_members, "high_confidence_member")
    add(args.candidate_representatives, "candidate_representative")
    if not records:
        raise ValueError("No protein sequences were prepared")

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    inventory = []
    with args.output_fasta.open("w") as handle:
        for query_id, record in records.items():
            handle.write(f">{query_id}\n")
            seq = record["sequence"]
            for start in range(0, len(seq), 60):
                handle.write(seq[start:start + 60] + "\n")
            inventory.append({
                "group": group,
                "query_id": query_id,
                "species": record["species"],
                "gene_id": record["gene_id"],
                "orthogroup": record["orthogroup"],
                "sequence_length_aa": len(seq),
                "sequence_md5": hashlib.md5(seq.encode()).hexdigest(),
                "is_high_confidence_member": "yes" if "high_confidence_member" in record["sources"] else "no",
                "is_candidate_representative": "yes" if "candidate_representative" in record["sources"] else "no",
            })
    write_tsv(args.inventory, list(inventory[0]), inventory)
    print(f"Prepared {len(records)} unique {group} proteins")


def read_emapper(path: Path):
    header = None
    rows = {}
    with path.open() as handle:
        for raw in handle:
            if raw.startswith("##"):
                continue
            if raw.startswith("#query"):
                header = raw.rstrip("\n").lstrip("#").split("\t")
                continue
            if raw.startswith("#") or not raw.strip():
                continue
            if header is None:
                raise ValueError(f"Could not locate #query header in {path}")
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < len(header):
                fields.extend([""] * (len(header) - len(fields)))
            row = dict(zip(header, fields))
            query = row.get("query", "")
            if not query:
                raise ValueError(f"Annotation row lacks query ID in {path}")
            if query in rows:
                raise ValueError(f"Duplicate eggNOG annotation row: {query}")
            rows[query] = row
    if header is None:
        raise ValueError(f"No eggNOG annotation header found in {path}")
    return header, rows


def first_existing(row, *names):
    for name in names:
        if name in row:
            return normal(row[name])
    return ""


def consensus(values):
    clean = [normal(v) for v in values if normal(v)]
    if not clean:
        return "", 0, 0.0, 0
    counts = Counter(clean)
    value, count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
    return value, count, count / len(clean), len(counts)


def union_terms(rows, field):
    terms = set()
    for row in rows:
        terms.update(split_terms(row.get(field, "")))
    return ",".join(sorted(terms))


def annotate_family_rows(group, ranked_rows, protein_rows_by_og):
    output = []
    summary_fields = [
        "n_queries", "n_seed_ortholog", "n_functionally_annotated", "annotation_fraction",
        "consensus_preferred_name", "preferred_name_support_n", "preferred_name_support_fraction",
        "unique_preferred_names", "consensus_description", "description_support_n",
        "description_support_fraction", "unique_descriptions", "cog_categories", "go_terms",
        "ec_numbers", "kegg_kos", "kegg_pathways", "kegg_modules", "kegg_reactions",
        "cazy_families", "bigg_reactions", "pfam_domains",
    ]
    for family in ranked_rows:
        orthogroup = family["orthogroup"]
        members = protein_rows_by_og.get(orthogroup, [])
        n = len(members)
        n_seed = sum(bool(normal(row.get("seed_ortholog", ""))) for row in members)
        functional_fields = ("Description", "Preferred_name", "GOs", "KEGG_ko", "PFAMs")
        n_annot = sum(any(normal(row.get(field, "")) for field in functional_fields) for row in members)
        pname, pname_n, pname_fraction, pname_unique = consensus(row.get("Preferred_name", "") for row in members)
        desc, desc_n, desc_fraction, desc_unique = consensus(row.get("Description", "") for row in members)
        extra = {
            "n_queries": n,
            "n_seed_ortholog": n_seed,
            "n_functionally_annotated": n_annot,
            "annotation_fraction": f"{(n_annot / n):.6f}" if n else "0.000000",
            "consensus_preferred_name": pname,
            "preferred_name_support_n": pname_n,
            "preferred_name_support_fraction": f"{pname_fraction:.6f}",
            "unique_preferred_names": pname_unique,
            "consensus_description": desc,
            "description_support_n": desc_n,
            "description_support_fraction": f"{desc_fraction:.6f}",
            "unique_descriptions": desc_unique,
            "cog_categories": union_terms(members, "COG_category"),
            "go_terms": union_terms(members, "GOs"),
            "ec_numbers": union_terms(members, "EC"),
            "kegg_kos": union_terms(members, "KEGG_ko"),
            "kegg_pathways": union_terms(members, "KEGG_Pathway"),
            "kegg_modules": union_terms(members, "KEGG_Module"),
            "kegg_reactions": union_terms(members, "KEGG_Reaction"),
            "cazy_families": union_terms(members, "CAZy"),
            "bigg_reactions": union_terms(members, "BiGG_Reaction"),
            "pfam_domains": union_terms(members, "PFAMs"),
        }
        merged = dict(family)
        merged.update(extra)
        output.append(merged)
    return output, summary_fields


def write_long_table(outdir, protein_rows, field, filename, term_label):
    rows = []
    seen = set()
    for row in protein_rows:
        for term in split_terms(row.get(field, "")):
            key = (row["query_id"], term)
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "group": row["group"], "orthogroup": row["orthogroup"],
                "species": row["species"], "gene_id": row["gene_id"],
                "query_id": row["query_id"], term_label: term,
            })
    fields = ["group", "orthogroup", "species", "gene_id", "query_id", term_label]
    write_tsv(outdir / filename, fields, rows)

    by_term = defaultdict(lambda: {"queries": set(), "orthogroups": set(), "species": set()})
    for row in rows:
        item = by_term[row[term_label]]
        item["queries"].add(row["query_id"])
        item["orthogroups"].add(row["orthogroup"])
        item["species"].add(row["species"])
    counts = [{
        "group": protein_rows[0]["group"] if protein_rows else "",
        term_label: term,
        "protein_count": len(value["queries"]),
        "orthogroup_count": len(value["orthogroups"]),
        "species_count": len(value["species"]),
    } for term, value in by_term.items()]
    counts.sort(key=lambda row: (-row["protein_count"], row[term_label]))
    count_fields = ["group", term_label, "protein_count", "orthogroup_count", "species_count"]
    write_tsv(outdir / filename.replace("_long.tsv", "_counts.tsv"), count_fields, counts)


def cmd_summarise(args):
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    inventory = read_tsv(args.inventory)
    inv_by_query = {row["query_id"]: row for row in inventory}
    if len(inv_by_query) != len(inventory):
        raise ValueError("Duplicate query IDs in inventory")
    emapper_header, emapper_by_query = read_emapper(args.annotations)
    unknown = sorted(set(emapper_by_query) - set(inv_by_query))
    if unknown:
        raise ValueError(f"eggNOG output contains {len(unknown)} unknown query IDs")

    canonical = [
        "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category",
        "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway",
        "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
        "BiGG_Reaction", "PFAMs",
    ]
    protein_rows = []
    for inv in inventory:
        query = inv["query_id"]
        ann = emapper_by_query.get(query, {})
        row = dict(inv)
        for field in canonical:
            row[field] = first_existing(ann, field)
        has_seed = bool(row["seed_ortholog"])
        has_function = any(row[field] for field in ("Description", "Preferred_name", "GOs", "KEGG_ko", "PFAMs"))
        row["annotation_status"] = "functional_annotation" if has_function else ("seed_ortholog_only" if has_seed else "unannotated")
        protein_rows.append(row)
    protein_fields = list(inventory[0]) + ["annotation_status"] + canonical
    write_tsv(outdir / "protein_functional_annotations.tsv", protein_fields, protein_rows)

    metric_fields = {
        "seed_ortholog": "seed_ortholog", "functional_annotation": None,
        "preferred_name": "Preferred_name", "description": "Description", "GO": "GOs",
        "EC": "EC", "KEGG_KO": "KEGG_ko", "KEGG_pathway": "KEGG_Pathway",
        "CAZy": "CAZy", "Pfam": "PFAMs",
    }
    coverage = []
    total = len(protein_rows)
    for metric, field in metric_fields.items():
        if field is None:
            count = sum(row["annotation_status"] == "functional_annotation" for row in protein_rows)
        else:
            count = sum(bool(normal(row.get(field, ""))) for row in protein_rows)
        coverage.append({
            "group": args.group, "metric": metric, "count": count,
            "total_queries": total, "fraction": f"{count / total:.6f}" if total else "0.000000",
        })
    write_tsv(outdir / "annotation_coverage_summary.tsv", list(coverage[0]), coverage)

    long_specs = [
        ("GOs", "go_annotations_long.tsv", "go_term"),
        ("EC", "ec_annotations_long.tsv", "ec_number"),
        ("KEGG_ko", "kegg_ko_annotations_long.tsv", "kegg_ko"),
        ("KEGG_Pathway", "kegg_pathway_annotations_long.tsv", "kegg_pathway"),
        ("CAZy", "cazy_annotations_long.tsv", "cazy_family"),
        ("PFAMs", "pfam_annotations_long.tsv", "pfam_domain"),
    ]
    for field, filename, label in long_specs:
        write_long_table(outdir, protein_rows, field, filename, label)

    ranked = read_tsv(args.ranked_families)
    proteins_by_og = defaultdict(list)
    for row in protein_rows:
        proteins_by_og[row["orthogroup"]].append(row)
    family_rows, summary_fields = annotate_family_rows(args.group, ranked, proteins_by_og)
    family_fields = list(ranked[0]) + summary_fields
    write_tsv(outdir / "candidate_family_functional_annotations.tsv", family_fields, family_rows)
    high_rows = [row for row in family_rows if row.get("candidate_tier") == "high"]
    write_tsv(outdir / "high_confidence_family_functional_annotations.tsv", family_fields, high_rows)

    family_by_og = {row["orthogroup"]: row for row in family_rows}
    terminal = read_tsv(args.terminal_events)
    terminal_out = []
    append_fields = [
        "n_queries", "n_functionally_annotated", "annotation_fraction", "consensus_preferred_name",
        "consensus_description", "go_terms", "kegg_kos", "kegg_pathways", "cazy_families", "pfam_domains",
    ]
    for event in terminal:
        merged = dict(event)
        family = family_by_og.get(event["orthogroup"], {})
        for field in append_fields:
            merged[field] = family.get(field, "")
        terminal_out.append(merged)
    write_tsv(outdir / "terminal_candidate_functional_annotations.tsv", list(terminal[0]) + append_fields, terminal_out)

    shutil.copy2(args.annotations, outdir / "raw_eggnog_annotations.tsv")
    if args.seed_orthologs and args.seed_orthologs.exists():
        shutil.copy2(args.seed_orthologs, outdir / "raw_eggnog_seed_orthologs.tsv")

    methods = f"""Step 7 Functional Annotation ({args.group})

Input design
- Every protein from Step 6 high-confidence families was included.
- Longest-per-species representatives from all CAFE-significant candidate families were added.
- Duplicate query identifiers were removed only when their amino-acid sequences were identical.

Annotation
- eggNOG-mapper 2.1.12, eggNOG database 5.0.2, DIAMOND mode.
- Protein input; eukaryota taxonomic scope; sensitive iterative DIAMOND search; DIAMOND algorithm selected automatically for compatibility.
- GO transfer restricted to non-electronic evidence; orthology-based Pfam transfer used without sequence realignment.
- The SQLite annotation database was read from disk; --dbmem was not used because jobs are capped at 16 GB RAM.

Interpretation limits
- Functional labels are orthology-based predictions, not experimental validation.
- Pfam labels in this run are transferred through eggNOG orthology and are not a substitute for direct InterProScan/HMMER domain validation.
- Family consensus is the most frequent non-missing annotation among queried family proteins; heterogeneous annotations are retained in protein-level and union-term columns.
- Bivalvia and Gastropoda orthogroup identifiers originate from separate OrthoFinder analyses and are not equivalent across groups.
"""
    (outdir / "METHODS_AND_LIMITATIONS.txt").write_text(methods)
    manifest = []
    for path in sorted(outdir.iterdir()):
        if path.is_file() and path.name not in {"COMPLETE", "file_manifest.tsv"}:
            manifest.append({"file": path.name, "bytes": path.stat().st_size})
    write_tsv(outdir / "file_manifest.tsv", ["file", "bytes"], manifest)
    (outdir / "COMPLETE").write_text("complete\n")
    print(f"Summarised {len(protein_rows)} proteins and {len(family_rows)} candidate families for {args.group}")


def cmd_combine(args):
    rows = []
    for group, group_dir in (("bivalvia", args.bivalvia), ("gastropoda", args.gastropoda)):
        if not (group_dir / "COMPLETE").exists():
            raise ValueError(f"Incomplete Step 7 output: {group_dir}")
        for row in read_tsv(group_dir / "annotation_coverage_summary.tsv"):
            rows.append(row)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "step07_annotation_coverage_summary.tsv", list(rows[0]), rows)
    run_rows = [
        {"group": "bivalvia", "result_directory": str(args.bivalvia.resolve())},
        {"group": "gastropoda", "result_directory": str(args.gastropoda.resolve())},
    ]
    write_tsv(args.output_dir / "step07_result_directories.tsv", list(run_rows[0]), run_rows)
    (args.output_dir / "COMPLETE").write_text("complete\n")
    print("Combined Bivalvia and Gastropoda Step 7 summaries")


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    prepare = sub.add_parser("prepare")
    prepare.add_argument("--group", required=True, choices=["bivalvia", "gastropoda"])
    prepare.add_argument("--high-confidence-members", type=Path, required=True)
    prepare.add_argument("--candidate-representatives", type=Path, required=True)
    prepare.add_argument("--output-fasta", type=Path, required=True)
    prepare.add_argument("--inventory", type=Path, required=True)
    prepare.set_defaults(func=cmd_prepare)

    summarise = sub.add_parser("summarise")
    summarise.add_argument("--group", required=True, choices=["bivalvia", "gastropoda"])
    summarise.add_argument("--inventory", type=Path, required=True)
    summarise.add_argument("--annotations", type=Path, required=True)
    summarise.add_argument("--seed-orthologs", type=Path)
    summarise.add_argument("--ranked-families", type=Path, required=True)
    summarise.add_argument("--terminal-events", type=Path, required=True)
    summarise.add_argument("--output-dir", type=Path, required=True)
    summarise.set_defaults(func=cmd_summarise)

    combine = sub.add_parser("combine")
    combine.add_argument("--bivalvia", type=Path, required=True)
    combine.add_argument("--gastropoda", type=Path, required=True)
    combine.add_argument("--output-dir", type=Path, required=True)
    combine.set_defaults(func=cmd_combine)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
