#!/usr/bin/env python3
"""Step 14: transparent integrated ranking of candidate families and genes."""

from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path("/path/to/mollusc_project")
BASE = ROOT / "results_or_reports/integrated_candidate_prioritisation/evidence_weighted_v1"
GROUPS = ("bivalvia", "gastropoda")

TREE = ROOT / "results_or_reports/candidate_family_phylogenetics/validated_top30_v1/latest/gene_tree_summary.tsv"
HYPHY = ROOT / "results_or_reports/molecular_evolution/hyphy_busted_absrel_relax_v1/latest/family_molecular_evolution_results.tsv"
EXPR_FAMILY = ROOT / "results_or_reports/transcriptomic_integration/single_library_salmon_v1/latest/orthogroup_transcriptomic_integration.tsv"
EXPR_GENE = ROOT / "results_or_reports/transcriptomic_integration/single_library_salmon_v1/latest/candidate_expression.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    with path.open(newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def f(value: str | None, default: float = 0.0) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def i(value: str | None, default: int = 0) -> int:
    try:
        return int(float(value)) if value not in (None, "") else default
    except ValueError:
        return default


def qsig(value: str | None) -> bool:
    return value not in (None, "") and f(value, 1.0) <= 0.05


def nonempty(row: dict[str, str], fields: tuple[str, ...]) -> bool:
    return any(row.get(field, "").strip() for field in fields)


def replace_link(link: Path, target: Path) -> None:
    if link.is_symlink():
        link.unlink()
    elif link.exists():
        raise RuntimeError(f"Refusing to replace non-symlink: {link}")
    link.symlink_to(target.name)


def annotation_path(group: str) -> Path:
    return ROOT / f"results_or_reports/functional_annotation/{group}/eggnog_mapper_v2_1_12/latest/candidate_family_functional_annotations.tsv"


def validation_family_path(group: str) -> Path:
    return ROOT / f"results_or_reports/candidate_sequence_validation/{group}/pfam_hmmer_v1/latest/family_validation_summary.tsv"


def validation_gene_path(group: str) -> Path:
    return ROOT / f"results_or_reports/candidate_sequence_validation/{group}/pfam_hmmer_v1/latest/sequence_validation.tsv"


def recon_path(group: str) -> Path:
    return ROOT / f"results_or_reports/gene_tree_species_tree_reconciliation/{group}/dl_lca_v1/latest/family_reconciliation_summary.tsv"


def enrichment_path(group: str) -> Path:
    return ROOT / f"results_or_reports/functional_enrichment/orthogroup_level_v1/latest/{group}/significant_enrichment_results.tsv"


def family_tier(score: float, breadth: int, risk: bool, validation_ok: float) -> str:
    if risk:
        return "manual_review_required"
    if score >= 75 and breadth >= 6 and validation_ok >= 0.50:
        return "tier1_high_confidence"
    if score >= 60 and breadth >= 5 and validation_ok >= 0.35:
        return "tier2_supported"
    if score >= 45 and breadth >= 4:
        return "tier3_provisional"
    return "lower_priority"


def gene_status_points(status: str) -> float:
    return {"pass": 20.0, "review": 14.0, "review_high": 6.0, "exclude_recommended": 0.0}.get(status, 5.0)


def expression_points(tpm: float) -> float:
    if tpm >= 10:
        return 10.0
    if tpm >= 1:
        return 8.0
    if tpm >= 0.1:
        return 5.0
    if tpm > 0:
        return 2.0
    return 0.0


def main(run_id: str) -> None:
    run = BASE / f"run_{run_id}"
    run.mkdir(parents=True, exist_ok=False)

    tree_rows = read_tsv(TREE)
    tree = {(r["group"], r["orthogroup"]): r for r in tree_rows}
    selected = set(tree)
    if len(selected) != 60:
        raise ValueError(f"Expected 60 selected group-family pairs, observed {len(selected)}")

    annotations = {}
    validations = {}
    recons = {}
    genes_by_group: dict[str, list[dict[str, str]]] = {}
    enrichment_rows: list[dict] = []
    for group in GROUPS:
        annotations.update({(group, r["orthogroup"]): r for r in read_tsv(annotation_path(group))})
        validations.update({(group, r["orthogroup"]): r for r in read_tsv(validation_family_path(group))})
        recons.update({(group, r["orthogroup"]): r for r in read_tsv(recon_path(group))})
        genes_by_group[group] = read_tsv(validation_gene_path(group))
        for r in read_tsv(enrichment_path(group)):
            enrichment_rows.append({
                "group": group,
                "comparison": r.get("comparison", ""),
                "target_definition": r.get("target_definition", ""),
                "ontology": r.get("ontology", ""),
                "namespace": r.get("namespace", ""),
                "term_id": r.get("term_id", ""),
                "term_name": r.get("term_name", ""),
                "fold_enrichment": r.get("fold_enrichment", ""),
                "q_bh": r.get("q_bh", ""),
                "interpretation": "context_only_not_scored",
            })

    hyphy = {(r["group"], r["orthogroup"]): r for r in read_tsv(HYPHY)}
    expr_family = {(r["group"], r["orthogroup"]): r for r in read_tsv(EXPR_FAMILY)}
    expr_gene = {(r["group"], r["species"], r["gene_id"]): r for r in read_tsv(EXPR_GENE)}

    family_rows: list[dict] = []
    family_by_key: dict[tuple[str, str], dict] = {}
    for key in sorted(selected):
        group, orthogroup = key
        tr = tree[key]
        ann = annotations.get(key, {})
        val = validations.get(key, {})
        rec = recons.get(key, {})
        mol = hyphy.get(key, {})
        exp = expr_family.get(key, {})

        rank = i(tr.get("rank"), 30)
        family_change = 25.0 * max(1.0 / 30.0, (31.0 - rank) / 30.0)

        ann_fraction = min(1.0, max(0.0, f(ann.get("annotation_fraction"))))
        functional = 4.0 * ann_fraction
        functional += 2.0 if ann.get("consensus_description", "").strip() else 0.0
        functional += 1.0 if ann.get("consensus_preferred_name", "").strip() else 0.0
        functional += 1.0 if ann.get("go_terms", "").strip() else 0.0
        functional += 1.0 if nonempty(ann, ("kegg_kos", "kegg_pathways", "ec_numbers")) else 0.0
        functional += 1.0 if ann.get("pfam_domains", "").strip() else 0.0

        size = max(1, i(val.get("family_size")))
        pass_n = i(val.get("pass_sequences"))
        review_n = i(val.get("review_sequences"))
        high_n = i(val.get("review_high_sequences"))
        exclude_n = i(val.get("exclude_recommended_sequences"))
        validation_ok = min(1.0, (pass_n + review_n) / size)
        pass_fraction = min(1.0, pass_n / size)
        direct_pfam_fraction = min(1.0, max(0.0, f(val.get("direct_pfam_fraction"))))
        high_risk_fraction = min(1.0, (high_n + exclude_n) / size)
        sequence_quality = 8.0 * validation_ok + 4.0 * pass_fraction + 3.0 * direct_pfam_fraction + 5.0 * (1.0 - high_risk_fraction)

        plastid = i(val.get("plastid_or_photosynthesis_flags"))
        te = i(val.get("transposable_element_flags"))
        viral = i(val.get("viral_or_phage_flags"))
        risk_fraction = (plastid + te + viral) / size
        critical_risk = exclude_n / size > 0.25 or risk_fraction > 0.50

        phylogenetic = 0.0
        if tr.get("status") == "complete":
            phylogenetic += 6.0
        if tr.get("support_mode") == "ufboot1000_alrt1000":
            phylogenetic += 4.0
        elif tr.get("support_mode", "").strip():
            phylogenetic += 2.0

        reconciliation = 0.0
        if rec:
            reconciliation += 4.0
            input_n = max(1, i(rec.get("input_sequences")))
            reconciliation += 3.0 * min(1.0, i(rec.get("reconciled_sequences")) / input_n)
            direction = tr.get("best_event_direction", "")
            if direction == "expansion" and i(rec.get("strong_support_duplications")) > 0:
                reconciliation += 3.0
            elif direction == "contraction" and i(rec.get("losses")) > 0:
                reconciliation += 3.0

        molecular_available = mol.get("execution_status") == "complete"
        molecular = 0.0
        if molecular_available:
            molecular += 8.0 if qsig(mol.get("busted_q_bh")) else 0.0
            molecular += 6.0 if i(mol.get("absrel_significant_branches_raw")) > 0 else 0.0
            molecular += 6.0 if qsig(mol.get("relax_q_bh")) else 0.0

        expr_available = bool(exp) and i(exp.get("species_with_expression_data")) > 0
        expression = 0.0
        expression_detection_fraction = 0.0
        if expr_available:
            copies = max(1, i(exp.get("candidate_copy_count")))
            expression_detection_fraction = min(1.0, i(exp.get("detected_copy_count_tpm_ge_0.1")) / copies)
            expression += 3.0 * expression_detection_fraction
            max_tpm = f(exp.get("max_within_library_tpm"))
            expression += 2.0 if max_tpm >= 10 else (1.5 if max_tpm >= 1 else (1.0 if max_tpm >= 0.1 else 0.0))

        raw = family_change + functional + sequence_quality + phylogenetic + reconciliation + molecular + expression
        available_max = 75.0 + (20.0 if molecular_available else 0.0) + (5.0 if expr_available else 0.0)
        normalized = 100.0 * raw / available_max
        domains = 4
        domains += int(functional > 0)
        domains += int(molecular_available)
        domains += int(expr_available)
        tier = family_tier(normalized, domains, critical_risk, validation_ok)

        evidence = {
            "group": group,
            "orthogroup": orthogroup,
            "group_rank_input": rank,
            "priority_tier": tier,
            "integrated_score_normalized": f"{normalized:.3f}",
            "raw_evidence_score": f"{raw:.3f}",
            "available_max_score": f"{available_max:.1f}",
            "evidence_domain_count": domains,
            "critical_risk_flag": "yes" if critical_risk else "no",
            "family_change_score_25": f"{family_change:.3f}",
            "functional_annotation_score_10": f"{functional:.3f}",
            "sequence_quality_score_20": f"{sequence_quality:.3f}",
            "gene_tree_score_10": f"{phylogenetic:.3f}",
            "reconciliation_score_10": f"{reconciliation:.3f}",
            "molecular_evolution_score_20": f"{molecular:.3f}" if molecular_available else "NA",
            "expression_support_score_5": f"{expression:.3f}" if expr_available else "NA",
            "best_event_node": tr.get("best_event_node", rec.get("cafe_best_event_node", "")),
            "best_event_direction": tr.get("best_event_direction", ""),
            "best_event_change": tr.get("best_event_change", ""),
            "original_selection_score": tr.get("selection_score", ""),
            "annotation_fraction": f"{ann_fraction:.6f}",
            "consensus_preferred_name": ann.get("consensus_preferred_name", ""),
            "consensus_description": ann.get("consensus_description", val.get("consensus_description", "")),
            "go_terms": ann.get("go_terms", ""),
            "kegg_kos": ann.get("kegg_kos", ""),
            "kegg_pathways": ann.get("kegg_pathways", ""),
            "pfam_domains": ann.get("pfam_domains", ""),
            "family_size": size,
            "validation_acceptable_fraction": f"{validation_ok:.6f}",
            "pass_fraction": f"{pass_fraction:.6f}",
            "direct_pfam_fraction": f"{direct_pfam_fraction:.6f}",
            "high_risk_fraction": f"{high_risk_fraction:.6f}",
            "exclude_recommended_sequences": exclude_n,
            "plastid_or_photosynthesis_flags": plastid,
            "transposable_element_flags": te,
            "viral_or_phage_flags": viral,
            "tree_status": tr.get("status", ""),
            "tree_support_mode": tr.get("support_mode", ""),
            "strong_support_duplications": rec.get("strong_support_duplications", ""),
            "reconciliation_losses": rec.get("losses", ""),
            "molecular_evolution_available": "yes" if molecular_available else "no",
            "busted_q_bh": mol.get("busted_q_bh", ""),
            "absrel_significant_branches_raw": mol.get("absrel_significant_branches_raw", ""),
            "relax_k": mol.get("relax_k", ""),
            "relax_q_bh": mol.get("relax_q_bh", ""),
            "relax_interpretation": mol.get("relax_interpretation", ""),
            "expression_data_available": "yes" if expr_available else "no",
            "species_with_expression_data": exp.get("species_with_expression_data", ""),
            "expression_detection_fraction": f"{expression_detection_fraction:.6f}" if expr_available else "",
            "max_within_library_tpm": exp.get("max_within_library_tpm", ""),
        }
        family_rows.append(evidence)
        family_by_key[key] = evidence

    # Rank within each group, then globally, with reproducible tie breaking.
    family_rows.sort(key=lambda r: (-f(r["integrated_score_normalized"]), r["group"], r["orthogroup"]))
    group_counter = defaultdict(int)
    for global_rank, row in enumerate(family_rows, 1):
        group_counter[row["group"]] += 1
        row["global_rank"] = global_rank
        row["group_rank_integrated"] = group_counter[row["group"]]

    family_fields = [
        "global_rank", "group", "group_rank_integrated", "orthogroup", "group_rank_input",
        "priority_tier", "integrated_score_normalized", "raw_evidence_score", "available_max_score",
        "evidence_domain_count", "critical_risk_flag", "family_change_score_25",
        "functional_annotation_score_10", "sequence_quality_score_20", "gene_tree_score_10",
        "reconciliation_score_10", "molecular_evolution_score_20", "expression_support_score_5",
        "best_event_node", "best_event_direction", "best_event_change", "original_selection_score",
        "annotation_fraction", "consensus_preferred_name", "consensus_description", "go_terms",
        "kegg_kos", "kegg_pathways", "pfam_domains", "family_size",
        "validation_acceptable_fraction", "pass_fraction", "direct_pfam_fraction",
        "high_risk_fraction", "exclude_recommended_sequences", "plastid_or_photosynthesis_flags",
        "transposable_element_flags", "viral_or_phage_flags", "tree_status", "tree_support_mode",
        "strong_support_duplications", "reconciliation_losses", "molecular_evolution_available",
        "busted_q_bh", "absrel_significant_branches_raw", "relax_k", "relax_q_bh",
        "relax_interpretation", "expression_data_available", "species_with_expression_data",
        "expression_detection_fraction", "max_within_library_tpm",
    ]
    write_tsv(run / "ranked_candidate_families.tsv", family_fields, family_rows)
    shortlisted = [r for r in family_rows if r["priority_tier"] in {"tier1_high_confidence", "tier2_supported"}]
    write_tsv(run / "high_confidence_family_shortlist.tsv", family_fields, shortlisted)

    # Gene ranking: family evidence 70%, sequence validation 20%, optional within-library expression 10%.
    gene_rows: list[dict] = []
    for group in GROUPS:
        for row in genes_by_group[group]:
            key = (group, row.get("orthogroup", ""))
            if key not in selected:
                continue
            family = family_by_key[key]
            family_component = 0.70 * f(family["integrated_score_normalized"])
            sequence_component = gene_status_points(row.get("validation_status", ""))
            exp = expr_gene.get((group, row.get("species", ""), row.get("gene_id", "")))
            exp_available = exp is not None
            exp_component = expression_points(f(exp.get("tpm"))) if exp else 0.0
            raw_gene = family_component + sequence_component + exp_component
            available_gene_max = 90.0 + (10.0 if exp_available else 0.0)
            gene_score = 100.0 * raw_gene / available_gene_max
            status = row.get("validation_status", "")
            if status == "exclude_recommended" or family["critical_risk_flag"] == "yes":
                tier = "manual_review_or_exclude"
            elif gene_score >= 75 and status == "pass":
                tier = "tier1_high_confidence"
            elif gene_score >= 60 and status in {"pass", "review"}:
                tier = "tier2_supported"
            elif gene_score >= 45 and status != "exclude_recommended":
                tier = "tier3_provisional"
            else:
                tier = "lower_priority"
            gene_rows.append({
                "group": group,
                "orthogroup": row.get("orthogroup", ""),
                "species": row.get("species", ""),
                "gene_id": row.get("gene_id", ""),
                "query_id": row.get("query_id", ""),
                "gene_priority_tier": tier,
                "gene_score_normalized": f"{gene_score:.3f}",
                "available_gene_max_score": f"{available_gene_max:.1f}",
                "family_integrated_score": family["integrated_score_normalized"],
                "family_priority_tier": family["priority_tier"],
                "sequence_validation_status": status,
                "sequence_validation_flags": row.get("validation_flags", ""),
                "sequence_length_aa": row.get("sequence_length_aa", ""),
                "direct_pfam_coverage_fraction": row.get("direct_pfam_coverage_fraction", ""),
                "direct_pfam_accessions": row.get("direct_pfam_accessions", ""),
                "direct_pfam_architecture": row.get("direct_pfam_architecture", ""),
                "description": row.get("description", ""),
                "preferred_name": row.get("preferred_name", ""),
                "expression_available": "yes" if exp_available else "no",
                "run_accession": exp.get("run_accession", "") if exp else "",
                "tpm": exp.get("tpm", "") if exp else "",
                "num_reads": exp.get("num_reads", "") if exp else "",
                "expression_support": exp.get("expression_support", "") if exp else "",
                "best_event_direction": family["best_event_direction"],
                "best_event_change": family["best_event_change"],
                "busted_q_bh": family["busted_q_bh"],
                "relax_q_bh": family["relax_q_bh"],
            })

    gene_rows.sort(key=lambda r: (-f(r["gene_score_normalized"]), r["group"], r["orthogroup"], r["species"], r["gene_id"]))
    gene_group_counter = defaultdict(int)
    for global_rank, row in enumerate(gene_rows, 1):
        gene_group_counter[row["group"]] += 1
        row["global_gene_rank"] = global_rank
        row["group_gene_rank"] = gene_group_counter[row["group"]]
    gene_fields = [
        "global_gene_rank", "group", "group_gene_rank", "orthogroup", "species", "gene_id",
        "query_id", "gene_priority_tier", "gene_score_normalized", "available_gene_max_score",
        "family_integrated_score", "family_priority_tier", "sequence_validation_status",
        "sequence_validation_flags", "sequence_length_aa", "direct_pfam_coverage_fraction",
        "direct_pfam_accessions", "direct_pfam_architecture", "description", "preferred_name",
        "expression_available", "run_accession", "tpm", "num_reads", "expression_support",
        "best_event_direction", "best_event_change", "busted_q_bh", "relax_q_bh",
    ]
    write_tsv(run / "ranked_candidate_genes.tsv", gene_fields, gene_rows)
    high_genes = [r for r in gene_rows if r["gene_priority_tier"] in {"tier1_high_confidence", "tier2_supported"}]
    write_tsv(run / "high_confidence_candidate_genes.tsv", gene_fields, high_genes)

    write_tsv(
        run / "significant_enrichment_context.tsv",
        ["group", "comparison", "target_definition", "ontology", "namespace", "term_id", "term_name", "fold_enrichment", "q_bh", "interpretation"],
        enrichment_rows,
    )

    summary_rows = []
    for group in GROUPS:
        fam = [r for r in family_rows if r["group"] == group]
        gen = [r for r in gene_rows if r["group"] == group]
        summary_rows.append({
            "group": group,
            "ranked_families": len(fam),
            "tier1_families": sum(r["priority_tier"] == "tier1_high_confidence" for r in fam),
            "tier2_families": sum(r["priority_tier"] == "tier2_supported" for r in fam),
            "manual_review_families": sum(r["priority_tier"] == "manual_review_required" for r in fam),
            "ranked_genes": len(gen),
            "tier1_genes": sum(r["gene_priority_tier"] == "tier1_high_confidence" for r in gen),
            "tier2_genes": sum(r["gene_priority_tier"] == "tier2_supported" for r in gen),
            "manual_review_or_exclude_genes": sum(r["gene_priority_tier"] == "manual_review_or_exclude" for r in gen),
            "families_with_molecular_evolution": sum(r["molecular_evolution_available"] == "yes" for r in fam),
            "families_with_expression_data": sum(r["expression_data_available"] == "yes" for r in fam),
        })
    write_tsv(
        run / "step14_summary.tsv",
        ["group", "ranked_families", "tier1_families", "tier2_families", "manual_review_families", "ranked_genes", "tier1_genes", "tier2_genes", "manual_review_or_exclude_genes", "families_with_molecular_evolution", "families_with_expression_data"],
        summary_rows,
    )

    rubric = """Step 14 integrated candidate prioritisation: scoring rubric

Family-level domains (maximum 100 when all evidence is available)
------------------------------------------------------------------
Family change evidence 25: rank-scaled among the preselected top 30 families per group.
Functional annotation 10: annotation coverage, description/name, GO, KEGG/EC and Pfam.
Sequence quality 20: pass/review fraction, pass fraction, direct-Pfam fraction and absence of
high-risk/exclude flags.
Gene-tree evidence 10: completed tree plus UFBoot1000/SH-aLRT1000 support mode.
Reconciliation 10: reconciliation availability, retained sequence fraction and a direction-
consistent supported duplication (expansion) or inferred loss (contraction).
Molecular evolution 20: BUSTED FDR (8), at least one aBSREL significant branch (6), RELAX FDR (6).
Transcriptomics 5: within-library candidate detection fraction (3) and maximum TPM support (2).

Missingness and interpretation
------------------------------
Scores are normalised by the maximum available evidence. Missing RNA-seq does not mean zero
expression, and families filtered from HyPhy for insufficient codon data are not assigned a
false negative molecular-evolution score. Evidence-domain count remains visible so equal scores
with different evidence breadth can be distinguished. Transcriptomic evidence is low weight
because each species has one selected library and is partly circular for transcript-derived
proteomes. Cross-species TPM values are not directly compared.

Risk gating and tiers
---------------------
Families with >25% exclude-recommended sequences or >50% combined plastid/photosynthesis,
transposable-element and viral/phage flags are assigned manual_review_required regardless of
score. Tier 1 requires score >=75, >=6 evidence domains and validation acceptable fraction >=0.5.
Tier 2 requires score >=60, >=5 domains and validation acceptable fraction >=0.35. Tier 3
requires score >=45 and >=4 domains. Gene ranking uses family evidence (70), sequence validation
(20) and optional within-library expression support (10), again normalised by available evidence.

Functional enrichment
---------------------
Significant Step 8 terms are exported as context but are not added to the score, because Step 8
used a candidate-internal background and adding the same annotations would double count evidence.

This ranking prioritises follow-up; it is not a probability of biological truth.
"""
    (run / "SCORING_RUBRIC.md").write_text(rubric)
    methods = """Step 14 integrates Steps 6-13 using deterministic, auditable rules. Inputs are
resolved from successful latest result links. No imputation of missing expression or molecular-
evolution results is performed. Rankings are generated at both orthogroup and sequence levels.
All source fields required to reproduce each score are retained in the output tables.
"""
    (run / "METHODS_AND_LIMITATIONS.txt").write_text(methods)

    sources = [
        ("step07_functional_annotation_bivalvia", annotation_path("bivalvia")),
        ("step07_functional_annotation_gastropoda", annotation_path("gastropoda")),
        ("step08_functional_enrichment_bivalvia", enrichment_path("bivalvia")),
        ("step08_functional_enrichment_gastropoda", enrichment_path("gastropoda")),
        ("step09_family_validation_bivalvia", validation_family_path("bivalvia")),
        ("step09_family_validation_gastropoda", validation_family_path("gastropoda")),
        ("step09_gene_validation_bivalvia", validation_gene_path("bivalvia")),
        ("step09_gene_validation_gastropoda", validation_gene_path("gastropoda")),
        ("step10_gene_trees", TREE),
        ("step11_reconciliation_bivalvia", recon_path("bivalvia")),
        ("step11_reconciliation_gastropoda", recon_path("gastropoda")),
        ("step12_molecular_evolution", HYPHY),
        ("step13_family_expression", EXPR_FAMILY),
        ("step13_gene_expression", EXPR_GENE),
    ]
    write_tsv(run / "evidence_source_manifest.tsv", ["evidence_source", "source_path"], [{"evidence_source": name, "source_path": path.resolve()} for name, path in sources])

    inventory = []
    for path in sorted(run.iterdir()):
        if path.is_file() and path.name not in {"COMPLETE", "result_manifest.tsv"}:
            inventory.append({"file": path.name, "bytes": path.stat().st_size})
    write_tsv(run / "result_manifest.tsv", ["file", "bytes"], inventory)
    (run / "COMPLETE").write_text(f"completed_run\t{run_id}\nranked_families\t{len(family_rows)}\nranked_genes\t{len(gene_rows)}\n")
    replace_link(BASE / "latest", run)
    print(f"Step 14 completed: {run}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 31_step14_integrated_prioritisation.py RUN_ID")
    main(sys.argv[1])
