#!/usr/bin/env python3
"""Prioritise CAFE5 Base-model candidate gene families for downstream annotation."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


GROUPS = ("bivalvia", "gastropoda")
NODE_RE = re.compile(r"^(?:(.+))?<([0-9]+)>$")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--family-fdr", type=float, default=0.05)
    parser.add_argument("--high-branch-p", type=float, default=0.01)
    parser.add_argument("--candidate-branch-p", type=float, default=0.05)
    parser.add_argument("--high-min-change", type=int, default=2)
    return parser.parse_args()


def tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        return list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def bh_adjust(records: dict[str, dict[str, object]]) -> None:
    ranked = sorted(records, key=lambda key: (float(records[key]["family_p"]), key))
    count = len(ranked)
    running = 1.0
    for reverse_rank, key in enumerate(reversed(ranked), start=1):
        rank = count - reverse_rank + 1
        raw = float(records[key]["family_p"]) * count / rank
        running = min(running, raw)
        records[key]["family_q_bh"] = min(1.0, running)


def parse_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def safe_neg_log10(value: float) -> float:
    return min(50.0, -math.log10(max(value, 1e-50)))


def read_fasta_selected(path: Path, wanted: set[str]) -> dict[str, str]:
    found: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current in wanted:
                    found[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if current in wanted:
            found[current] = "".join(chunks)
    return found


def write_fasta(path: Path, entries: list[tuple[str, str, str, str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for group, family, species, gene, sequence in entries:
            handle.write(f">{species}|{gene}|{family} group={group}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")


def load_matrix(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or reader.fieldnames[0] != "FamilyID":
            raise ValueError(f"Unexpected CAFE matrix header: {path}")
        columns = reader.fieldnames[1:]
        rows = {row["FamilyID"]: row for row in reader}
    return columns, rows


def analyse_group(
    project: Path,
    output: Path,
    group: str,
    family_fdr: float,
    high_branch_p: float,
    candidate_branch_p: float,
    high_min_change: int,
) -> dict[str, object]:
    cafe = project / "results_or_reports" / "gene_family_evolution" / group / "relative_time_screen_v1" / "latest"
    ortho = project / "orthofinder_results" / "by_class" / f"{group}_v3_1_5" / "Results_Aug04" / "Orthogroups"
    distribution_path = project / "results_or_reports" / "orthogroup_distribution" / "latest" / group / "orthogroup_distribution.tsv"
    proteomes = project / "orthofinder_primary_proteomes"
    required = [
        cafe / "cafe_base_error" / "Base_family_results.txt",
        cafe / "cafe_base_error" / "Base_change.tab",
        cafe / "cafe_base_error" / "Base_count.tab",
        cafe / "cafe_base_error" / "Base_branch_probabilities.tab",
        ortho / "Orthogroups.GeneCount.tsv",
        ortho / "Orthogroups.tsv",
        distribution_path,
    ]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    group_out = output / group
    group_out.mkdir()

    families: dict[str, dict[str, object]] = {}
    for row in tsv_rows(cafe / "cafe_base_error" / "Base_family_results.txt"):
        family = row["#FamilyID"]
        families[family] = {
            "group": group,
            "orthogroup": family,
            "family_p": float(row["pvalue"]),
            "cafe_significant_0.05": row["Significant at 0.05"].lower(),
        }
    bh_adjust(families)

    node_columns, changes = load_matrix(cafe / "cafe_base_error" / "Base_change.tab")
    count_columns, reconstructed = load_matrix(cafe / "cafe_base_error" / "Base_count.tab")
    probability_columns, probabilities = load_matrix(cafe / "cafe_base_error" / "Base_branch_probabilities.tab")
    if node_columns != count_columns or node_columns != probability_columns:
        raise ValueError(f"CAFE node columns differ for {group}")
    if set(families) != set(changes) or set(families) != set(reconstructed):
        raise ValueError(f"CAFE family identifiers differ between family/change/count tables for {group}")
    # CAFE5 intentionally writes branch probabilities only for family-wide
    # significant families, not for every family in the likelihood analysis.
    significant_ids = {
        family for family, values in families.items()
        if str(values["cafe_significant_0.05"]) == "y"
    }
    missing_probabilities = significant_ids - set(probabilities)
    if missing_probabilities:
        preview = ",".join(sorted(missing_probabilities)[:10])
        raise ValueError(f"Missing branch probabilities for {len(missing_probabilities)} significant families in {group}: {preview}")

    distribution = {row["orthogroup"]: row for row in tsv_rows(distribution_path)}
    gene_counts = {row["Orthogroup"]: row for row in tsv_rows(ortho / "Orthogroups.GeneCount.tsv")}
    if not set(families).issubset(distribution) or not set(families).issubset(gene_counts):
        raise ValueError(f"Missing orthogroup metadata for {group}")
    species = [column for column in next(iter(gene_counts.values())) if column not in {"Orthogroup", "Total"}]

    all_events: list[dict[str, object]] = []
    event_by_family: dict[str, list[dict[str, object]]] = defaultdict(list)
    for family in sorted(families):
        family_p = float(families[family]["family_p"])
        family_q = float(families[family]["family_q_bh"])
        if str(families[family]["cafe_significant_0.05"]) != "y":
            continue
        for label in node_columns:
            change = int(changes[family][label])
            branch_p = parse_float(probabilities[family][label])
            if change == 0 or branch_p is None:
                continue
            child_count = int(reconstructed[family][label])
            parent_count = child_count - change
            match = NODE_RE.match(label)
            if match is None:
                raise ValueError(f"Cannot parse CAFE node label: {label}")
            terminal_species = match.group(1) or ""
            is_terminal = bool(terminal_species)
            observed_count: int | str = ""
            if is_terminal:
                if terminal_species not in species:
                    raise ValueError(f"Unknown terminal species {terminal_species} in {group}")
                observed_count = int(gene_counts[family][terminal_species])
            direction = "expansion" if change > 0 else "contraction"
            fold_magnitude = ((child_count + 1) / (parent_count + 1)) if change > 0 else ((parent_count + 1) / (child_count + 1))
            if family_q <= family_fdr and branch_p <= high_branch_p and abs(change) >= high_min_change:
                tier = "high"
            elif branch_p <= candidate_branch_p:
                tier = "candidate"
            else:
                tier = "supporting"
            score = (
                safe_neg_log10(family_q)
                + safe_neg_log10(branch_p)
                + math.log2(1 + abs(change))
                + math.log2(max(1.0, fold_magnitude))
                + (1.0 if is_terminal else 0.0)
            )
            event = {
                "group": group,
                "orthogroup": family,
                "family_p": f"{family_p:.12g}",
                "family_q_bh": f"{family_q:.12g}",
                "node": label,
                "node_id": match.group(2),
                "terminal_species": terminal_species,
                "is_terminal": "yes" if is_terminal else "no",
                "direction": direction,
                "change": change,
                "absolute_change": abs(change),
                "parent_reconstructed_count": parent_count,
                "child_reconstructed_count": child_count,
                "observed_terminal_count": observed_count,
                "fold_magnitude_pseudocount1": f"{fold_magnitude:.6f}",
                "branch_probability": f"{branch_p:.12g}",
                "event_tier": tier,
                "evidence_score": f"{score:.6f}",
            }
            all_events.append(event)
            event_by_family[family].append(event)

    tier_order = {"high": 0, "candidate": 1, "supporting": 2}
    all_events.sort(key=lambda row: (tier_order[str(row["event_tier"])], -float(row["evidence_score"]), str(row["orthogroup"]), str(row["node"])))
    event_columns = [
        "group", "orthogroup", "family_p", "family_q_bh", "node", "node_id", "terminal_species", "is_terminal",
        "direction", "change", "absolute_change", "parent_reconstructed_count", "child_reconstructed_count",
        "observed_terminal_count", "fold_magnitude_pseudocount1", "branch_probability", "event_tier", "evidence_score",
    ]
    write_tsv(group_out / "significant_family_branch_changes.tsv", all_events, event_columns)
    high_events = [row for row in all_events if row["event_tier"] == "high"]
    terminal_events = [row for row in all_events if row["is_terminal"] == "yes" and row["event_tier"] in {"high", "candidate"}]
    write_tsv(group_out / "high_confidence_branch_events.tsv", high_events, event_columns)
    write_tsv(group_out / "terminal_candidate_events.tsv", terminal_events, event_columns)

    family_rows: list[dict[str, object]] = []
    high_families: set[str] = set()
    raw_candidate_families: set[str] = set()
    for family in sorted(families):
        meta = families[family]
        family_p = float(meta["family_p"])
        family_q = float(meta["family_q_bh"])
        if str(meta["cafe_significant_0.05"]) != "y":
            continue
        raw_candidate_families.add(family)
        events = event_by_family.get(family, [])
        high = [event for event in events if event["event_tier"] == "high"]
        candidates = [event for event in events if event["event_tier"] == "candidate"]
        if high:
            tier = "high"
            high_families.add(family)
        elif family_q <= family_fdr and candidates:
            tier = "fdr_supported"
        elif candidates:
            tier = "candidate"
        else:
            tier = "family_only"
        best = max(events, key=lambda event: float(event["evidence_score"])) if events else None
        dist = distribution[family]
        counts = gene_counts[family]
        row: dict[str, object] = {
            "group": group,
            "orthogroup": family,
            "family_p": f"{family_p:.12g}",
            "family_q_bh": f"{family_q:.12g}",
            "candidate_tier": tier,
            "distribution_category": dist["category"],
            "species_present": dist["n_species_present"],
            "occupancy_fraction": dist["occupancy_fraction"],
            "total_observed_copies": counts["Total"],
            "max_observed_copies": dist["max_copies"],
            "max_copy_species": dist["max_copy_species"],
            "high_events": len(high),
            "candidate_events": len(candidates),
            "terminal_expansion_events": sum(1 for event in events if event["is_terminal"] == "yes" and event["direction"] == "expansion" and event["event_tier"] in {"high", "candidate"}),
            "terminal_contraction_events": sum(1 for event in events if event["is_terminal"] == "yes" and event["direction"] == "contraction" and event["event_tier"] in {"high", "candidate"}),
            "best_event_node": best["node"] if best else "",
            "best_event_direction": best["direction"] if best else "",
            "best_event_change": best["change"] if best else "",
            "best_event_branch_probability": best["branch_probability"] if best else "",
            "evidence_score": best["evidence_score"] if best else f"{safe_neg_log10(family_q):.6f}",
            "functional_annotation_status": "pending",
        }
        for name in species:
            row[name] = counts[name]
        family_rows.append(row)
    family_tier_order = {"high": 0, "fdr_supported": 1, "candidate": 2, "family_only": 3}
    family_rows.sort(key=lambda row: (family_tier_order[str(row["candidate_tier"])], -float(row["evidence_score"]), str(row["orthogroup"])))
    family_columns = [
        "group", "orthogroup", "family_p", "family_q_bh", "candidate_tier", "distribution_category", "species_present",
        "occupancy_fraction", "total_observed_copies", "max_observed_copies", "max_copy_species", "high_events",
        "candidate_events", "terminal_expansion_events", "terminal_contraction_events", "best_event_node",
        "best_event_direction", "best_event_change", "best_event_branch_probability", "evidence_score",
        "functional_annotation_status", *species,
    ]
    write_tsv(group_out / "candidate_families_ranked.tsv", family_rows, family_columns)
    write_tsv(group_out / "high_confidence_families.tsv", [row for row in family_rows if row["candidate_tier"] == "high"], family_columns)

    membership_rows = tsv_rows(ortho / "Orthogroups.tsv")
    member_map: dict[str, dict[str, list[str]]] = {}
    membership_output: list[dict[str, object]] = []
    for row in membership_rows:
        family = row["Orthogroup"]
        if family not in raw_candidate_families:
            continue
        member_map[family] = {}
        tier = next(item["candidate_tier"] for item in family_rows if item["orthogroup"] == family)
        for name in species:
            genes = [gene.strip() for gene in row[name].split(",") if gene.strip()]
            member_map[family][name] = genes
            for gene in genes:
                membership_output.append({"group": group, "orthogroup": family, "candidate_tier": tier, "species": name, "gene_id": gene})
    write_tsv(group_out / "candidate_gene_membership.tsv", membership_output, ["group", "orthogroup", "candidate_tier", "species", "gene_id"])

    high_terminal_expansions = {
        (str(event["orthogroup"]), str(event["terminal_species"]))
        for event in high_events
        if event["is_terminal"] == "yes" and event["direction"] == "expansion"
    }
    needed_by_species: dict[str, set[str]] = defaultdict(set)
    terminal_needed: set[tuple[str, str, str]] = set()
    for family in high_families:
        for name, genes in member_map.get(family, {}).items():
            needed_by_species[name].update(genes)
            if (family, name) in high_terminal_expansions:
                terminal_needed.update((family, name, gene) for gene in genes)

    sequences: dict[tuple[str, str], str] = {}
    missing: list[dict[str, object]] = []
    for name in species:
        fasta = proteomes / f"{name}.faa"
        if not fasta.is_file() or fasta.stat().st_size == 0:
            raise FileNotFoundError(fasta)
        found = read_fasta_selected(fasta, needed_by_species[name])
        for gene in sorted(needed_by_species[name]):
            if gene in found:
                sequences[(name, gene)] = found[gene]
            else:
                missing.append({"group": group, "species": name, "gene_id": gene})
    write_tsv(group_out / "missing_high_confidence_sequences.tsv", missing, ["group", "species", "gene_id"])
    if missing:
        raise RuntimeError(f"Missing {len(missing)} high-confidence protein sequences for {group}")

    all_entries: list[tuple[str, str, str, str, str]] = []
    terminal_entries: list[tuple[str, str, str, str, str]] = []
    representative_entries: list[tuple[str, str, str, str, str]] = []
    sequence_rows: list[dict[str, object]] = []
    for family in sorted(high_families):
        for name in species:
            genes = member_map.get(family, {}).get(name, [])
            present: list[tuple[str, str]] = []
            for gene in genes:
                sequence = sequences[(name, gene)]
                present.append((gene, sequence))
                all_entries.append((group, family, name, gene, sequence))
                if (family, name, gene) in terminal_needed:
                    terminal_entries.append((group, family, name, gene, sequence))
                sequence_rows.append({
                    "group": group,
                    "orthogroup": family,
                    "species": name,
                    "gene_id": gene,
                    "protein_length_aa": len(sequence.rstrip("*")),
                    "terminal_high_confidence_expansion_member": "yes" if (family, name, gene) in terminal_needed else "no",
                })
            if present:
                gene, sequence = max(present, key=lambda item: (len(item[1].rstrip("*")), item[0]))
                representative_entries.append((group, family, name, gene, sequence))
    write_fasta(group_out / "high_confidence_family_all_members.faa", all_entries)
    write_fasta(group_out / "high_confidence_terminal_expansion_proteins.faa", terminal_entries)
    write_fasta(group_out / "annotation_representatives_longest_per_species.faa", representative_entries)
    write_tsv(
        group_out / "high_confidence_sequence_inventory.tsv",
        sequence_rows,
        ["group", "orthogroup", "species", "gene_id", "protein_length_aa", "terminal_high_confidence_expansion_member"],
    )

    summary = {
        "group": group,
        "cafe_base_families": len(families),
        "cafe_significant_families": sum(1 for value in families.values() if str(value["cafe_significant_0.05"]) == "y"),
        "bh_fdr_le_0.05_families": sum(1 for value in families.values() if float(value["family_q_bh"]) <= family_fdr),
        "high_confidence_families": len(high_families),
        "high_confidence_branch_events": len(high_events),
        "high_confidence_terminal_expansions": sum(1 for event in high_events if event["is_terminal"] == "yes" and event["direction"] == "expansion"),
        "high_confidence_terminal_contractions": sum(1 for event in high_events if event["is_terminal"] == "yes" and event["direction"] == "contraction"),
        "candidate_gene_members": len(membership_output),
        "high_confidence_family_proteins": len(all_entries),
        "terminal_expansion_proteins": len(terminal_entries),
        "annotation_representatives": len(representative_entries),
        "missing_high_confidence_sequences": len(missing),
    }
    write_tsv(group_out / "summary.tsv", [summary], list(summary))
    return summary


def main() -> None:
    args = arguments()
    project = args.project.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    summaries = [
        analyse_group(
            project,
            output,
            group,
            args.family_fdr,
            args.high_branch_p,
            args.candidate_branch_p,
            args.high_min_change,
        )
        for group in GROUPS
    ]
    write_tsv(output / "step06_summary.tsv", summaries, list(summaries[0]))
    methods = [
        "Step 6: Candidate gene-family identification",
        "CAFE5 Base-model family-wide p-values were used for both groups.",
        "Benjamini-Hochberg correction was applied independently within each class-specific CAFE analysis.",
        f"High-confidence branch event: family FDR <= {args.family_fdr}, branch probability <= {args.high_branch_p}, and absolute reconstructed change >= {args.high_min_change}.",
        f"Candidate branch event: CAFE5 family significant flag = y, branch probability <= {args.candidate_branch_p}, and non-zero reconstructed change.",
        "CAFE terminal reconstructed counts, observed OrthoFinder counts, orthogroup distribution, and protein membership were retained as separate evidence fields.",
        "Gastropoda Gamma results were excluded because all input families reported >20% failure rates.",
        "Bivalvia and Gastropoda orthogroup identifiers come from separate OrthoFinder analyses and must not be matched across groups by ID alone.",
        "Functional annotation is pending; these tables identify statistical/copy-number candidates, not biological function or adaptation.",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
    ]
    (output / "METHODS_AND_LIMITATIONS.txt").write_text("\n".join(methods) + "\n", encoding="utf-8")
    (output / "COMPLETE").write_text(datetime.now(timezone.utc).isoformat() + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
