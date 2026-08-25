#!/usr/bin/env python3
"""Step 8: orthogroup-level functional over-representation analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
from collections import defaultdict
from pathlib import Path


MISSING = {"", "-", "NA", "N/A", "None", "none", "null"}
ONTOLOGY_FIELDS = {
    "GO": "go_terms",
    "KEGG_KO": "kegg_kos",
    "KEGG_PATHWAY": "kegg_pathways",
    "PFAM": "pfam_domains",
    "CAZY": "cazy_families",
    "EC": "ec_numbers",
}


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_terms(value):
    value = (value or "").strip()
    if value in MISSING:
        return []
    return sorted({term.strip() for term in value.split(",") if term.strip() and term.strip() not in MISSING})


def parse_go_obo(path: Path | None):
    terms = {}
    alt_to_primary = {}
    if path is None or not path.exists():
        return terms, alt_to_primary
    current = None
    records = []
    with path.open(errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                if current:
                    records.append(current)
                current = {"alt_id": []}
            elif line.startswith("["):
                if current:
                    records.append(current)
                current = None
            elif current is not None and ": " in line:
                key, value = line.split(": ", 1)
                if key == "alt_id":
                    current["alt_id"].append(value)
                elif key in {"id", "name", "namespace", "is_obsolete"}:
                    current[key] = value
    if current:
        records.append(current)
    for record in records:
        go_id = record.get("id")
        if not go_id or record.get("is_obsolete") == "true":
            continue
        terms[go_id] = {
            "name": record.get("name", ""),
            "namespace": record.get("namespace", "unknown"),
        }
        for alt_id in record.get("alt_id", []):
            alt_to_primary[alt_id] = go_id
    return terms, alt_to_primary


def log_choose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_upper(k, population, successes, draws):
    upper = min(successes, draws)
    if k > upper:
        return 0.0
    lower = max(0, draws - (population - successes))
    start = max(k, lower)
    logs = [
        log_choose(successes, x)
        + log_choose(population - successes, draws - x)
        - log_choose(population, draws)
        for x in range(start, upper + 1)
    ]
    peak = max(logs)
    return min(1.0, math.exp(peak) * sum(math.exp(value - peak) for value in logs))


def bh_adjust(rows):
    if not rows:
        return
    ordered = sorted(enumerate(rows), key=lambda item: item[1]["p_value_num"])
    m = len(rows)
    adjusted = [1.0] * m
    running = 1.0
    for rank_index in range(m - 1, -1, -1):
        original_index, row = ordered[rank_index]
        rank = rank_index + 1
        running = min(running, row["p_value_num"] * m / rank)
        adjusted[original_index] = min(1.0, running)
    for row, value in zip(rows, adjusted):
        row["q_bh_num"] = value


def safe_label(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_").lower()


def build_comparisons(family_rows, terminal_rows, min_species_set):
    all_families = {row["orthogroup"] for row in family_rows}
    comparisons = {
        "high_confidence": {
            "definition": "Step 6 high-confidence families versus all CAFE-significant candidate families",
            "families": {row["orthogroup"] for row in family_rows if row.get("candidate_tier") == "high"},
        },
    }
    high_terminal = [row for row in terminal_rows if row.get("event_tier") == "high"]
    for direction in ("expansion", "contraction"):
        comparisons[f"terminal_{direction}_high"] = {
            "definition": f"Families with a high-confidence terminal {direction} versus all CAFE-significant candidate families",
            "families": {row["orthogroup"] for row in high_terminal if row.get("direction") == direction},
        }
    species_events = defaultdict(set)
    for row in high_terminal:
        species = row.get("terminal_species", "")
        direction = row.get("direction", "")
        if species and direction in {"expansion", "contraction"}:
            species_events[(species, direction)].add(row["orthogroup"])
    for (species, direction), families in sorted(species_events.items()):
        if len(families) >= min_species_set:
            key = f"species_{safe_label(species)}_{direction}"
            comparisons[key] = {
                "definition": f"High-confidence terminal {direction} families in {species} versus all CAFE-significant candidate families",
                "families": families,
            }
    for value in comparisons.values():
        value["families"] &= all_families
    return all_families, comparisons


def analyse_group(group, group_dir, go_terms, go_alt, args, run_dir):
    families_path = group_dir / "candidate_family_functional_annotations.tsv"
    terminal_path = group_dir / "terminal_candidate_functional_annotations.tsv"
    if not (group_dir / "COMPLETE").exists():
        raise ValueError(f"Step 7 output is incomplete: {group_dir}")
    family_rows = read_tsv(families_path)
    terminal_rows = read_tsv(terminal_path)
    family_by_og = {row["orthogroup"]: row for row in family_rows}
    if len(family_by_og) != len(family_rows):
        raise ValueError(f"Duplicate orthogroups in {families_path}")
    all_families, comparisons = build_comparisons(family_rows, terminal_rows, args.min_species_set)

    candidate_set_rows = []
    comparison_meta = {}
    for name, spec in comparisons.items():
        comparison_meta[name] = {
            "group": group,
            "comparison": name,
            "target_definition": spec["definition"],
            "target_family_count_all": len(spec["families"]),
            "background_family_count_all": len(all_families),
        }
        for orthogroup in sorted(spec["families"]):
            candidate_set_rows.append({
                "group": group, "comparison": name,
                "target_definition": spec["definition"], "orthogroup": orthogroup,
            })

    all_results = []
    comparison_summaries = []
    for ontology, field in ONTOLOGY_FIELDS.items():
        term_to_families = defaultdict(set)
        term_details = {}
        families_with_terms = set()
        for orthogroup, row in family_by_og.items():
            raw_terms = split_terms(row.get(field, ""))
            canonical_terms = []
            for term in raw_terms:
                canonical = go_alt.get(term, term) if ontology == "GO" else term
                canonical_terms.append(canonical)
                if ontology == "GO":
                    detail = go_terms.get(canonical, {})
                    term_details[canonical] = (detail.get("name", ""), detail.get("namespace", "unknown"))
                else:
                    term_details.setdefault(canonical, ("", ""))
            if canonical_terms:
                families_with_terms.add(orthogroup)
            for term in set(canonical_terms):
                term_to_families[term].add(orthogroup)
        universe = all_families & families_with_terms
        population = len(universe)

        for comparison, spec in comparisons.items():
            target_all = spec["families"]
            target = target_all & universe
            draws = len(target)
            rows_for_test = []
            if draws >= args.min_annotated_target and population > draws:
                for term, members in term_to_families.items():
                    background_members = members & universe
                    target_members = background_members & target
                    successes = len(background_members)
                    observed = len(target_members)
                    if successes < args.min_background_count or observed < args.min_target_count:
                        continue
                    p_value = hypergeom_upper(observed, population, successes, draws)
                    target_fraction = observed / draws
                    background_fraction = successes / population
                    fold = target_fraction / background_fraction if background_fraction else float("inf")
                    a = observed
                    b = draws - observed
                    c = successes - observed
                    d = population - successes - b
                    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
                    term_name, namespace = term_details.get(term, ("", ""))
                    rows_for_test.append({
                        "group": group,
                        "comparison": comparison,
                        "target_definition": spec["definition"],
                        "ontology": ontology,
                        "namespace": namespace,
                        "term_id": term,
                        "term_name": term_name,
                        "background_annotated_families": population,
                        "target_families_all": len(target_all),
                        "target_annotated_families": draws,
                        "background_term_families": successes,
                        "target_term_families": observed,
                        "background_fraction": background_fraction,
                        "target_fraction": target_fraction,
                        "fold_enrichment": fold,
                        "odds_ratio_haldane": odds,
                        "p_value_num": p_value,
                    })
            strata = defaultdict(list)
            for row in rows_for_test:
                stratum = row["namespace"] if ontology == "GO" else ontology
                strata[stratum].append(row)
            for stratum_rows in strata.values():
                bh_adjust(stratum_rows)
            all_results.extend(rows_for_test)
            significant = sum(row.get("q_bh_num", 1.0) <= 0.05 for row in rows_for_test)
            comparison_summaries.append({
                **comparison_meta[comparison],
                "ontology": ontology,
                "background_annotated_families": population,
                "target_annotated_families": draws,
                "tested_terms": len(rows_for_test),
                "significant_terms_fdr05": significant,
            })

    all_results.sort(key=lambda row: (row["comparison"], row["ontology"], row["namespace"], row.get("q_bh_num", 1), row["p_value_num"], row["term_id"]))
    result_fields = [
        "group", "comparison", "target_definition", "ontology", "namespace", "term_id", "term_name",
        "background_annotated_families", "target_families_all", "target_annotated_families",
        "background_term_families", "target_term_families", "background_fraction", "target_fraction",
        "fold_enrichment", "odds_ratio_haldane", "p_value", "q_bh", "significant_fdr05",
    ]
    formatted = []
    for row in all_results:
        out = {key: row.get(key, "") for key in result_fields}
        out["background_fraction"] = f"{row['background_fraction']:.8g}"
        out["target_fraction"] = f"{row['target_fraction']:.8g}"
        out["fold_enrichment"] = f"{row['fold_enrichment']:.8g}"
        out["odds_ratio_haldane"] = f"{row['odds_ratio_haldane']:.8g}"
        out["p_value"] = f"{row['p_value_num']:.12g}"
        out["q_bh"] = f"{row.get('q_bh_num', 1.0):.12g}"
        out["significant_fdr05"] = "yes" if row.get("q_bh_num", 1.0) <= 0.05 else "no"
        formatted.append(out)

    group_out = run_dir / group
    write_tsv(group_out / "all_enrichment_results.tsv", result_fields, formatted)
    write_tsv(group_out / "significant_enrichment_results.tsv", result_fields, [row for row in formatted if row["significant_fdr05"] == "yes"])
    summary_fields = list(comparison_summaries[0]) if comparison_summaries else ["group", "comparison"]
    write_tsv(group_out / "comparison_summary.tsv", summary_fields, comparison_summaries)
    set_fields = list(candidate_set_rows[0]) if candidate_set_rows else ["group", "comparison", "target_definition", "orthogroup"]
    write_tsv(group_out / "candidate_sets.tsv", set_fields, candidate_set_rows)

    for ontology in ONTOLOGY_FIELDS:
        subset = [row for row in formatted if row["ontology"] == ontology]
        write_tsv(group_out / f"{ontology.lower()}_enrichment.tsv", result_fields, subset)
    for namespace in ("biological_process", "molecular_function", "cellular_component", "unknown"):
        subset = [row for row in formatted if row["ontology"] == "GO" and row["namespace"] == namespace]
        write_tsv(group_out / f"go_{namespace}_enrichment.tsv", result_fields, subset)

    return {
        "group": group,
        "background_candidate_families": len(all_families),
        "comparisons": len(comparisons),
        "tested_term_rows": len(formatted),
        "significant_term_rows_fdr05": sum(row["significant_fdr05"] == "yes" for row in formatted),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bivalvia", type=Path, required=True)
    parser.add_argument("--gastropoda", type=Path, required=True)
    parser.add_argument("--go-obo", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-target-count", type=int, default=3)
    parser.add_argument("--min-background-count", type=int, default=5)
    parser.add_argument("--min-annotated-target", type=int, default=5)
    parser.add_argument("--min-species-set", type=int, default=5)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    go_terms, go_alt = parse_go_obo(args.go_obo)
    summaries = []
    summaries.append(analyse_group("bivalvia", args.bivalvia, go_terms, go_alt, args, args.output_dir))
    summaries.append(analyse_group("gastropoda", args.gastropoda, go_terms, go_alt, args, args.output_dir))
    write_tsv(args.output_dir / "step08_summary.tsv", list(summaries[0]), summaries)

    go_checksum = ""
    if args.go_obo and args.go_obo.exists():
        go_checksum = hashlib.sha256(args.go_obo.read_bytes()).hexdigest()
    methods = f"""Step 8 Functional Enrichment Analysis

Statistical unit and background
- The orthogroup/family, not the individual protein copy, is the statistical unit. This prevents expanded families from receiving extra weight merely because they contain more paralogues.
- Analyses are performed separately for Bivalvia and Gastropoda because their OrthoFinder orthogroup identifiers are independent.
- The tested background is all CAFE-significant candidate families from Step 6 that received at least one term in the relevant Step 7 ontology.
- This conditional background supports prioritisation within the candidate set; it is not a whole-proteome enrichment test.

Candidate sets
- Step 6 high-confidence families.
- Families with high-confidence terminal expansions or contractions.
- Species-specific high-confidence terminal expansion/contraction sets containing at least {args.min_species_set} families.

Terms and testing
- GO, KEGG Orthology, KEGG pathways, Pfam, CAZy and EC terms transferred by eggNOG-mapper in Step 7.
- One-sided hypergeometric over-representation test.
- Terms require at least {args.min_target_count} target families and {args.min_background_count} background families.
- Comparisons require at least {args.min_annotated_target} annotated target families.
- Benjamini-Hochberg correction is applied separately within each group, comparison and ontology; GO namespaces are corrected separately.
- Significance threshold: BH FDR <= 0.05.

GO ontology
- GO metadata source: {args.go_obo if args.go_obo else 'not supplied'}
- Parsed current non-obsolete GO terms: {len(go_terms)}
- SHA256: {go_checksum}

Limitations
- Enrichment reflects orthology-transferred functional predictions and requires later biological interpretation and sequence/domain validation.
- Because only candidate families were annotated in Step 7, these results identify functions over-represented among the strongest candidates relative to other CAFE-significant candidates, not relative to every gene in each proteome.
- Very small lineage sets are omitted, and significant terms should not be interpreted as independent because ontology terms overlap hierarchically.
"""
    (args.output_dir / "METHODS_AND_LIMITATIONS.txt").write_text(methods)
    manifest = []
    for path in sorted(args.output_dir.rglob("*")):
        if path.is_file() and path.name not in {"COMPLETE", "file_manifest.tsv"}:
            manifest.append({"file": str(path.relative_to(args.output_dir)), "bytes": path.stat().st_size})
    write_tsv(args.output_dir / "file_manifest.tsv", ["file", "bytes"], manifest)
    (args.output_dir / "COMPLETE").write_text("complete\n")
    print("Step 8 functional enrichment complete")


if __name__ == "__main__":
    main()
