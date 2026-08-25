#!/usr/bin/env python3
"""Step 15: sensitivity and bias analysis for integrated candidate rankings.

This is deliberately a low-file-count analysis because the BluePebble work
fileset is close to its inode soft quota.  It operates on the 60 family rows
from Step 14 and writes compact, audit-friendly TSV reports.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


COMPONENTS = {
    "change": ("family_change_score_25", 25.0),
    "annotation": ("functional_annotation_score_10", 10.0),
    "sequence": ("sequence_quality_score_20", 20.0),
    "tree": ("gene_tree_score_10", 10.0),
    "reconciliation": ("reconciliation_score_10", 10.0),
    "evolution": ("molecular_evolution_score_20", 20.0),
    "expression": ("expression_support_score_5", 5.0),
}

SCENARIOS = [
    ("baseline", "Original Step 14 evidence-weighted score and risk gate."),
    ("no_rna", "Remove the low-weight transcriptomic component; missing RNA remains non-penalising."),
    ("no_hyphy", "Remove molecular-evolution evidence from all families."),
    ("no_function_annotation", "Remove functional-annotation evidence from all families."),
    ("phylogeny_sequence_only", "Use sequence validation, gene-tree and reconciliation evidence only."),
    ("strict_sequence_quality", "Replace the sequence component with a stricter pass/Pfam/risk formulation."),
    ("strict_contamination_risk", "Place any family with exclude, plastid, TE or viral flags below clean families."),
    ("family_change_dominant", "Increase family-change weight to 40% and down-weight annotation/structure."),
    ("annotation_dominant", "Increase annotation weight to 30% to test annotation-source leverage."),
    ("complete_case_only", "Place families lacking any of seven evidence domains below complete-case families."),
    ("eggnog_annotation_only", "Derive annotation contribution only from eggNOG family annotation coverage."),
    ("direct_pfam_annotation_only", "Derive annotation contribution only from direct Pfam family coverage."),
    ("proteome_quality_adjusted", "Down-weight family-change evidence by BUSCO completeness and duplication."),
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def num(value, default: float | None = 0.0) -> float | None:
    if value is None or str(value).strip() in {"", "NA", "NaN", "nan"}:
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def yes(value) -> bool:
    return str(value).strip().lower() in {"yes", "true", "1", "complete", "completed"}


def component_ratio(row: dict[str, str], name: str) -> float | None:
    field, maximum = COMPONENTS[name]
    available = True
    if name == "evolution":
        available = yes(row.get("molecular_evolution_available"))
    elif name == "expression":
        available = yes(row.get("expression_data_available"))
    if not available:
        return None
    return max(0.0, min(1.0, num(row.get(field), 0.0) / maximum))


def weighted_score(row: dict[str, str], weights: dict[str, float], overrides=None) -> float:
    overrides = overrides or {}
    raw = 0.0
    maximum = 0.0
    for name, weight in weights.items():
        ratio = overrides.get(name, component_ratio(row, name))
        if ratio is None:
            continue
        raw += weight * max(0.0, min(1.0, ratio))
        maximum += weight
    return 100.0 * raw / maximum if maximum else 0.0


def risk_counts(row: dict[str, str]) -> tuple[int, int, int, int]:
    return tuple(int(num(row.get(field), 0.0)) for field in (
        "exclude_recommended_sequences",
        "plastid_or_photosynthesis_flags",
        "transposable_element_flags",
        "viral_or_phage_flags",
    ))


def match_event_species(node: str, species: dict[str, dict[str, str]]) -> str:
    node = node or ""
    matches = [sp for sp in species if node == sp or node.startswith(sp + "<")]
    return max(matches, key=len) if matches else ""


def quality_factor(row: dict[str, str], species: dict[str, dict[str, str]], group_means: dict[str, tuple[float, float]]) -> tuple[float, str, float, float]:
    event_species = match_event_species(row.get("best_event_node", ""), species)
    if event_species:
        complete = num(species[event_species].get("complete_pct"), 0.0)
        duplicated = num(species[event_species].get("duplicated_pct"), 0.0)
    else:
        complete, duplicated = group_means[row["group"]]
    factor = max(0.0, min(1.0, complete / 100.0)) * max(0.0, 1.0 - duplicated / 100.0)
    return factor, event_species, complete, duplicated


def scenario_result(row: dict[str, str], scenario: str, species, group_means) -> tuple[float, bool, str]:
    standard = {name: maximum for name, (_, maximum) in COMPONENTS.items()}
    base = num(row.get("integrated_score_normalized"), 0.0)
    critical = str(row.get("critical_risk_flag", "no")).lower() == "yes"
    if scenario == "baseline":
        return base, not critical, "step14_score"
    if scenario == "no_rna":
        return weighted_score(row, {k: v for k, v in standard.items() if k != "expression"}), not critical, "expression_removed"
    if scenario == "no_hyphy":
        return weighted_score(row, {k: v for k, v in standard.items() if k != "evolution"}), not critical, "molecular_evolution_removed"
    if scenario == "no_function_annotation":
        return weighted_score(row, {k: v for k, v in standard.items() if k != "annotation"}), not critical, "functional_annotation_removed"
    if scenario == "phylogeny_sequence_only":
        weights = {"sequence": 20.0, "tree": 10.0, "reconciliation": 10.0}
        return weighted_score(row, weights), not critical, "sequence_tree_reconciliation_only"
    if scenario == "strict_sequence_quality":
        acceptable = num(row.get("validation_acceptable_fraction"), 0.0)
        passed = num(row.get("pass_fraction"), 0.0)
        pfam = num(row.get("direct_pfam_fraction"), 0.0)
        high_risk = num(row.get("high_risk_fraction"), 0.0)
        family_size = max(1.0, num(row.get("family_size"), 1.0))
        excluded = num(row.get("exclude_recommended_sequences"), 0.0) / family_size
        strict_ratio = max(0.0, min(1.0, 0.50 * acceptable + 0.30 * passed + 0.20 * pfam - 0.50 * high_risk - 0.50 * excluded))
        eligible = (not critical and acceptable >= 0.75 and excluded == 0 and high_risk <= 0.10)
        return weighted_score(row, standard, {"sequence": strict_ratio}), eligible, "strict_sequence_formula_and_gate"
    if scenario == "strict_contamination_risk":
        clean = sum(risk_counts(row)) == 0
        return base, clean and not critical, "zero_high_risk_flags_required"
    if scenario == "family_change_dominant":
        weights = {"change": 40, "annotation": 5, "sequence": 15, "tree": 8, "reconciliation": 8, "evolution": 20, "expression": 4}
        return weighted_score(row, weights), not critical, "change_weight_40"
    if scenario == "annotation_dominant":
        weights = {"change": 15, "annotation": 30, "sequence": 15, "tree": 8, "reconciliation": 8, "evolution": 20, "expression": 4}
        return weighted_score(row, weights), not critical, "annotation_weight_30"
    if scenario == "complete_case_only":
        complete = int(num(row.get("evidence_domain_count"), 0.0)) >= 7
        return base, complete and not critical, "seven_domains_required"
    if scenario == "eggnog_annotation_only":
        annotation_ratio = num(row.get("annotation_fraction"), 0.0)
        return weighted_score(row, standard, {"annotation": annotation_ratio}), not critical, "annotation_from_eggnog_coverage"
    if scenario == "direct_pfam_annotation_only":
        annotation_ratio = num(row.get("direct_pfam_fraction"), 0.0)
        return weighted_score(row, standard, {"annotation": annotation_ratio}), not critical, "annotation_from_direct_pfam_coverage"
    if scenario == "proteome_quality_adjusted":
        factor, event_species, complete, duplicated = quality_factor(row, species, group_means)
        change = component_ratio(row, "change")
        note = f"event_species={event_species or 'group_mean'};C={complete:.1f};D={duplicated:.1f};factor={factor:.4f}"
        return weighted_score(row, standard, {"change": (change or 0.0) * factor}), not critical, note
    raise ValueError(scenario)


def rank_values(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        rank = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = rank
        i = j + 1
    return ranks


def pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2:
        return float("nan")
    mx, my = statistics.mean(xs), statistics.mean(ys)
    numerator = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return numerator / (dx * dy) if dx and dy else float("nan")


def spearman(xs: list[float], ys: list[float]) -> float:
    return pearson(rank_values(xs), rank_values(ys))


def tier(score: float, domains: int, eligible: bool) -> str:
    if not eligible:
        return "manual_review_or_filtered"
    if score >= 75 and domains >= 6:
        return "tier1_high_confidence"
    if score >= 60 and domains >= 5:
        return "tier2_supported"
    if score >= 45 and domains >= 4:
        return "tier3_exploratory"
    return "lower_priority"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--step14-run", required=True, type=Path)
    parser.add_argument("--busco", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)

    families = read_tsv(args.step14_run / "ranked_candidate_families.tsv")
    busco_rows = read_tsv(args.busco)
    species = {row["species"]: row for row in busco_rows}
    group_species = {
        "bivalvia": {"Acanthocardia_echinata", "Americardia_media", "Cerastoderma_edule", "Fragum_sueziense", "Tridacna_gigas", "Tridacna_maxima"},
        "gastropoda": {"Aplysia_californica", "Berghia_stephanieae", "Bullacta_exarata", "Elysia_chlorotica", "Elysia_crispata", "Elysia_marginata", "Onchidella_celtica", "Plakobranchus_ocellatus"},
    }
    group_means = {}
    for group, members in group_species.items():
        selected = [species[s] for s in members if s in species]
        group_means[group] = (
            statistics.mean(num(r.get("complete_pct"), 0.0) for r in selected),
            statistics.mean(num(r.get("duplicated_pct"), 0.0) for r in selected),
        )

    rankings = []
    for scenario, _ in SCENARIOS:
        for group in sorted({row["group"] for row in families}):
            subset = []
            for row in families:
                if row["group"] != group:
                    continue
                score, eligible, note = scenario_result(row, scenario, species, group_means)
                subset.append((row, score, eligible, note))
            # The baseline must reproduce Step 14 exactly.  Alternative strict
            # gates may place ineligible families below eligible families, but
            # reapplying that gate to baseline would move Step 14 manual-review
            # families and make the baseline-versus-itself correlation < 1.
            if scenario == "baseline":
                subset.sort(key=lambda item: int(num(item[0].get("group_rank_integrated"), 0.0)))
            else:
                subset.sort(key=lambda item: (not item[2], -item[1], item[0]["orthogroup"]))
            for rank, (row, score, eligible, note) in enumerate(subset, 1):
                baseline_rank = int(num(row.get("group_rank_integrated"), 0.0))
                domains = int(num(row.get("evidence_domain_count"), 0.0))
                rankings.append({
                    "scenario": scenario,
                    "group": group,
                    "orthogroup": row["orthogroup"],
                    "scenario_rank": rank,
                    "baseline_rank": baseline_rank,
                    "rank_shift_vs_baseline": baseline_rank - rank,
                    "scenario_score": f"{score:.6f}",
                    "baseline_score": row.get("integrated_score_normalized", ""),
                    "scenario_eligible": "yes" if eligible else "no",
                    "scenario_tier": row.get("priority_tier", "") if scenario == "baseline" else tier(score, domains, eligible),
                    "baseline_tier": row.get("priority_tier", ""),
                    "evidence_domain_count": domains,
                    "critical_risk_flag": row.get("critical_risk_flag", ""),
                    "scenario_note": note,
                })

    ranking_fields = list(rankings[0].keys())
    write_tsv(args.out / "scenario_family_rankings.tsv", rankings, ranking_fields)

    by_key = {(r["scenario"], r["group"], r["orthogroup"]): r for r in rankings}
    comparison = []
    groups = sorted({r["group"] for r in families})
    for scenario, _ in SCENARIOS:
        for group in groups:
            current = [r for r in rankings if r["scenario"] == scenario and r["group"] == group]
            baseline_ranks = [int(r["baseline_rank"]) for r in current]
            scenario_ranks = [int(r["scenario_rank"]) for r in current]
            base10 = {r["orthogroup"] for r in current if int(r["baseline_rank"]) <= 10}
            base20 = {r["orthogroup"] for r in current if int(r["baseline_rank"]) <= 20}
            if scenario == "baseline":
                now10 = {r["orthogroup"] for r in current if int(r["scenario_rank"]) <= 10}
                now20 = {r["orthogroup"] for r in current if int(r["scenario_rank"]) <= 20}
            else:
                now10 = {r["orthogroup"] for r in current if int(r["scenario_rank"]) <= 10 and r["scenario_eligible"] == "yes"}
                now20 = {r["orthogroup"] for r in current if int(r["scenario_rank"]) <= 20 and r["scenario_eligible"] == "yes"}
            base_tier1 = {r["orthogroup"] for r in current if r["baseline_tier"] == "tier1_high_confidence"}
            now_tier1 = {r["orthogroup"] for r in current if r["scenario_tier"] == "tier1_high_confidence"}
            comparison.append({
                "scenario": scenario,
                "group": group,
                "families": len(current),
                "eligible_families": sum(r["scenario_eligible"] == "yes" for r in current),
                "spearman_rho_vs_baseline": f"{spearman(baseline_ranks, scenario_ranks):.6f}",
                "top10_overlap": len(base10 & now10),
                "top10_jaccard": f"{len(base10 & now10) / max(1, len(base10 | now10)):.6f}",
                "top20_overlap": len(base20 & now20),
                "top20_jaccard": f"{len(base20 & now20) / max(1, len(base20 | now20)):.6f}",
                "baseline_tier1": len(base_tier1),
                "scenario_tier1": len(now_tier1),
                "tier1_retained": len(base_tier1 & now_tier1),
                "maximum_absolute_rank_shift": max(abs(int(r["rank_shift_vs_baseline"])) for r in current),
            })
    write_tsv(args.out / "scenario_comparison_summary.tsv", comparison, list(comparison[0].keys()))

    stability = []
    for row in families:
        family_scenarios = [by_key[(scenario, row["group"], row["orthogroup"])] for scenario, _ in SCENARIOS]
        ranks = [int(r["scenario_rank"]) for r in family_scenarios]
        eligible_n = sum(r["scenario_eligible"] == "yes" for r in family_scenarios)
        top10_n = sum(int(r["scenario_rank"]) <= 10 and r["scenario_eligible"] == "yes" for r in family_scenarios)
        top20_n = sum(int(r["scenario_rank"]) <= 20 and r["scenario_eligible"] == "yes" for r in family_scenarios)
        max_shift = max(abs(int(r["rank_shift_vs_baseline"])) for r in family_scenarios)
        stability.append({
            "group": row["group"],
            "orthogroup": row["orthogroup"],
            "baseline_rank": row.get("group_rank_integrated", ""),
            "baseline_score": row.get("integrated_score_normalized", ""),
            "baseline_tier": row.get("priority_tier", ""),
            "mean_scenario_rank": f"{statistics.mean(ranks):.3f}",
            "rank_sd": f"{statistics.pstdev(ranks):.3f}",
            "best_rank": min(ranks),
            "worst_rank": max(ranks),
            "maximum_absolute_rank_shift": max_shift,
            "eligible_scenarios": eligible_n,
            "top10_scenarios": top10_n,
            "top20_scenarios": top20_n,
            "stability_class": "stable_high" if top10_n >= math.ceil(0.8 * len(SCENARIOS)) and max_shift <= 5 else ("sensitive" if max_shift >= 10 or eligible_n < len(SCENARIOS) - 2 else "moderately_stable"),
        })
    stability.sort(key=lambda r: (r["group"], float(r["mean_scenario_rank"]), r["orthogroup"]))
    write_tsv(args.out / "candidate_rank_stability.tsv", stability, list(stability[0].keys()))
    outliers = [r for r in stability if r["stability_class"] == "sensitive"]
    write_tsv(args.out / "rank_shift_outliers.tsv", outliers, list(stability[0].keys()))

    quality = []
    for row in families:
        factor, event_species, complete, duplicated = quality_factor(row, species, group_means)
        adjusted = by_key[("proteome_quality_adjusted", row["group"], row["orthogroup"])]
        quality.append({
            "group": row["group"],
            "orthogroup": row["orthogroup"],
            "best_event_node": row.get("best_event_node", ""),
            "matched_event_species": event_species or "group_mean_fallback",
            "busco_complete_pct": f"{complete:.3f}",
            "busco_duplicated_pct": f"{duplicated:.3f}",
            "quality_factor": f"{factor:.6f}",
            "baseline_rank": row.get("group_rank_integrated", ""),
            "quality_adjusted_rank": adjusted["scenario_rank"],
            "rank_shift_vs_baseline": adjusted["rank_shift_vs_baseline"],
            "baseline_score": row.get("integrated_score_normalized", ""),
            "quality_adjusted_score": adjusted["scenario_score"],
        })
    write_tsv(args.out / "proteome_quality_bias_assessment.tsv", quality, list(quality[0].keys()))

    thresholds = []
    for threshold in (80, 85, 90, 95):
        for group in groups:
            qrows = [r for r in quality if r["group"] == group]
            matched = [r for r in qrows if r["matched_event_species"] != "group_mean_fallback"]
            below = [r for r in matched if float(r["busco_complete_pct"]) < threshold]
            top10_below = [r for r in below if int(r["baseline_rank"]) <= 10]
            thresholds.append({
                "group": group,
                "busco_complete_threshold": threshold,
                "candidate_families": len(qrows),
                "event_species_matched": len(matched),
                "matched_families_below_threshold": len(below),
                "baseline_top10_below_threshold": len(top10_below),
                "below_threshold_orthogroups": ",".join(r["orthogroup"] for r in below) or "none",
            })
    write_tsv(args.out / "quality_threshold_sensitivity.tsv", thresholds, list(thresholds[0].keys()))

    annotation_bias = []
    for row in families:
        egg = by_key[("eggnog_annotation_only", row["group"], row["orthogroup"])]
        pfam = by_key[("direct_pfam_annotation_only", row["group"], row["orthogroup"])]
        annotation_bias.append({
            "group": row["group"],
            "orthogroup": row["orthogroup"],
            "baseline_rank": row.get("group_rank_integrated", ""),
            "eggnog_rank": egg["scenario_rank"],
            "direct_pfam_rank": pfam["scenario_rank"],
            "eggnog_minus_pfam_rank": int(egg["scenario_rank"]) - int(pfam["scenario_rank"]),
            "annotation_fraction": row.get("annotation_fraction", ""),
            "direct_pfam_fraction": row.get("direct_pfam_fraction", ""),
            "consensus_preferred_name": row.get("consensus_preferred_name", ""),
        })
    annotation_bias.sort(key=lambda r: (r["group"], -abs(int(r["eggnog_minus_pfam_rank"])), r["orthogroup"]))
    write_tsv(args.out / "annotation_source_bias.tsv", annotation_bias, list(annotation_bias[0].keys()))

    summary = []
    for group in groups:
        group_stability = [r for r in stability if r["group"] == group]
        group_comp = [r for r in comparison if r["group"] == group and r["scenario"] != "baseline"]
        summary.append({
            "group": group,
            "candidate_families": len(group_stability),
            "scenarios": len(SCENARIOS),
            "stable_high_families": sum(r["stability_class"] == "stable_high" for r in group_stability),
            "moderately_stable_families": sum(r["stability_class"] == "moderately_stable" for r in group_stability),
            "sensitive_families": sum(r["stability_class"] == "sensitive" for r in group_stability),
            "minimum_spearman_rho": f"{min(float(r['spearman_rho_vs_baseline']) for r in group_comp):.6f}",
            "minimum_top10_overlap": min(int(r["top10_overlap"]) for r in group_comp),
            "maximum_rank_shift": max(int(r["maximum_absolute_rank_shift"]) for r in group_comp),
        })
    write_tsv(args.out / "step15_summary.tsv", summary, list(summary[0].keys()))

    write_tsv(args.out / "SCENARIO_DEFINITIONS.tsv", [
        {"scenario": name, "definition": definition} for name, definition in SCENARIOS
    ], ["scenario", "definition"])

    methods = """Step 15 sensitivity and bias analysis

The 60 Step 14 family candidates were re-ranked under deterministic alternative evidence
weights, evidence-removal scenarios, stricter validation/risk gates, annotation-source
alternatives and BUSCO-based quality adjustment. Scores are normalised over evidence available
within each scenario. Missing transcriptomic or HyPhy evidence is never treated as a biological
negative. Scenario agreement is reported using Spearman rank correlation, top-10/top-20 overlap,
Jaccard similarity, tier-1 retention and maximum rank shift.

Proteome-quality analysis adjusts only the family-change component. It uses the BUSCO score of
the species that labels the best Step 6 event when this can be matched; otherwise it uses the
class mean. This is a diagnostic sensitivity analysis, not a re-run of OrthoFinder or CAFE with
different species sets. Likewise, eggNOG-only and direct-Pfam-only scenarios quantify annotation
source leverage and do not imply that either source is a gold standard.

Strict gates place affected families below eligible families but retain them in output for audit.
Rank stability is evidence robustness, not probability of biological truth. Follow-up should
prioritise families that remain highly ranked across scenarios and separately inspect sensitive
families for annotation, contamination or proteome-quality artefacts.
"""
    (args.out / "METHODS_AND_LIMITATIONS.txt").write_text(methods, encoding="utf-8")

    manifest_rows = []
    for path in sorted(args.out.iterdir()):
        if path.is_file():
            manifest_rows.append({"file": path.name, "bytes": path.stat().st_size})
    write_tsv(args.out / "result_manifest.tsv", manifest_rows, ["file", "bytes"])
    (args.out / "COMPLETE").write_text("Step 15 sensitivity and bias analysis complete\n", encoding="utf-8")
    print(f"Step 15 completed: {args.out}")


if __name__ == "__main__":
    main()
