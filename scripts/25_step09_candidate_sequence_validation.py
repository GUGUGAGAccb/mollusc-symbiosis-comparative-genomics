#!/usr/bin/env python3
"""Step 9 candidate-sequence validation and evidence-based triage."""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median


VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")
PLASTID_GO = {"GO:0009507", "GO:0009536", "GO:0009526", "GO:0015979", "GO:0019684"}
RISK_PATTERNS = {
    "plastid_or_photosynthesis": re.compile(
        r"chloroplast|plastid|photosystem|chlorophyll|light[- ]harvest|rubisco|ribulose[- ]?bisphosphate carboxylase|phycobil",
        re.I,
    ),
    "transposable_element": re.compile(
        r"reverse transcript|retrotrans|transposase|transposable element|retroviral|gag[- ]?pol|\bRVT[_-]|\bDDE[_-]?Tnp|integrase core",
        re.I,
    ),
    "viral_or_phage": re.compile(r"viral|virus|virion|capsid|bacteriophage|\bphage\b", re.I),
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
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before header in {path}")
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def write_fasta(path: Path, records, query_ids):
    with path.open("w") as handle:
        for query_id in query_ids:
            sequence = records[query_id]
            handle.write(f">{query_id}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")


def parse_domtbl(path: Path):
    hits = defaultdict(list)
    rows = []
    with path.open(errors="replace") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split(maxsplit=22)
            if len(fields) < 22:
                raise ValueError(f"Malformed HMMER domtblout row: {raw[:120]}")
            target_name, target_accession = fields[0], fields[1]
            query_id = fields[3]
            accession = target_accession.split(".")[0] if target_accession != "-" else target_name
            row = {
                "query_id": query_id,
                "pfam_name": target_name,
                "pfam_accession": accession,
                "full_evalue": fields[6],
                "full_score": fields[7],
                "domain_index": fields[9],
                "domain_total": fields[10],
                "conditional_evalue": fields[11],
                "independent_evalue": fields[12],
                "domain_score": fields[13],
                "hmm_from": int(fields[15]),
                "hmm_to": int(fields[16]),
                "alignment_from": int(fields[17]),
                "alignment_to": int(fields[18]),
                "envelope_from": int(fields[19]),
                "envelope_to": int(fields[20]),
                "accuracy": fields[21],
                "description": fields[22] if len(fields) > 22 else "",
            }
            hits[query_id].append(row)
            rows.append(row)
    for query_id in hits:
        hits[query_id].sort(key=lambda row: (row["alignment_from"], row["alignment_to"], row["pfam_accession"]))
    return hits, rows


def interval_coverage(intervals, length):
    if not intervals or length <= 0:
        return 0, 0.0
    merged = []
    for start, end in sorted(intervals):
        start = max(1, start)
        end = min(length, end)
        if end < start:
            continue
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    covered = sum(end - start + 1 for start, end in merged)
    return covered, covered / length


def sequence_entropy(sequence):
    residues = [aa for aa in sequence if aa in VALID_AA and aa != "*"]
    if not residues:
        return 0.0
    counts = Counter(residues)
    total = len(residues)
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def taxonomy_map(db_path: Path, taxids):
    result = {}
    if not db_path.exists():
        return result
    con = sqlite3.connect(str(db_path))
    merged = dict(con.execute("SELECT taxid_old, taxid_new FROM merged"))
    for original in sorted(set(taxids)):
        taxid = merged.get(original, original)
        row = con.execute("SELECT spname, rank, track FROM species WHERE taxid=?", (taxid,)).fetchone()
        if not row:
            result[original] = {"taxid": taxid, "species": "", "rank": "", "track": "", "class": "unknown"}
            continue
        species, rank, track = row
        lineage = {int(value) for value in str(track).split(",") if value.isdigit()}
        if 10239 in lineage:
            label = "virus"
        elif 2 in lineage:
            label = "bacteria"
        elif 2157 in lineage:
            label = "archaea"
        elif 33208 in lineage:
            label = "metazoa"
        elif 4751 in lineage:
            label = "fungi"
        elif 33090 in lineage:
            label = "viridiplantae"
        elif 33630 in lineage:
            label = "alveolata"
        elif 33634 in lineage:
            label = "stramenopiles"
        elif 2759 in lineage:
            label = "other_eukaryote"
        else:
            label = "unknown"
        result[original] = {"taxid": taxid, "species": species, "rank": rank, "track": track, "class": label}
    con.close()
    return result


def parse_seed_taxid(seed):
    match = re.match(r"^(\d+)[.]", (seed or "").strip())
    return int(match.group(1)) if match else None


def flag_text(row, direct_names):
    text = " ".join([
        row.get("Description", ""), row.get("Preferred_name", ""), row.get("PFAMs", ""),
        " ".join(direct_names),
    ])
    flags = {name for name, pattern in RISK_PATTERNS.items() if pattern.search(text)}
    if PLASTID_GO & set((row.get("GOs", "") or "").split(",")):
        flags.add("plastid_or_photosynthesis")
    return sorted(flags)


def cmd_validate(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    proteins = read_tsv(args.annotations)
    sequences = read_fasta(args.fasta)
    protein_by_query = {row["query_id"]: row for row in proteins}
    if set(protein_by_query) != set(sequences):
        missing_fasta = set(protein_by_query) - set(sequences)
        missing_table = set(sequences) - set(protein_by_query)
        raise ValueError(f"FASTA/table mismatch: missing_fasta={len(missing_fasta)}, missing_table={len(missing_table)}")
    pfam_hits, pfam_rows = parse_domtbl(args.domtblout)
    unknown_queries = set(pfam_hits) - set(sequences)
    if unknown_queries:
        raise ValueError(f"Pfam output has {len(unknown_queries)} unknown query IDs")

    family_input = read_tsv(args.family_table)
    family_meta = {row["orthogroup"]: row for row in family_input}
    family_queries = defaultdict(list)
    md5_queries = defaultdict(list)
    taxids = []
    for row in proteins:
        family_queries[row["orthogroup"]].append(row["query_id"])
        md5_queries[row["sequence_md5"]].append(row["query_id"])
        taxid = parse_seed_taxid(row.get("seed_ortholog", ""))
        if taxid is not None:
            taxids.append(taxid)
    taxa = taxonomy_map(args.taxonomy_db, taxids)

    family_stats = {}
    for orthogroup, query_ids in family_queries.items():
        lengths = [len(sequences[query]) for query in query_ids]
        with_domains = [query for query in query_ids if pfam_hits.get(query)]
        domain_counts = Counter()
        for query in with_domains:
            domain_counts.update({hit["pfam_accession"] for hit in pfam_hits[query]})
        core_domains = set()
        if len(with_domains) >= 3:
            threshold = max(2, math.ceil(0.60 * len(with_domains)))
            core_domains = {domain for domain, count in domain_counts.items() if count >= threshold}
        family_stats[orthogroup] = {
            "median_length": median(lengths),
            "n": len(query_ids),
            "n_with_domains": len(with_domains),
            "domain_rich_fraction": len(with_domains) / len(query_ids),
            "core_domains": core_domains,
        }

    duplicate_group_rows = []
    duplicate_context = {}
    duplicate_index = 0
    for md5, query_ids in sorted(md5_queries.items()):
        if len(query_ids) < 2:
            continue
        duplicate_index += 1
        orthogroups = {protein_by_query[q]["orthogroup"] for q in query_ids}
        species = {protein_by_query[q]["species"] for q in query_ids}
        context = "cross_orthogroup" if len(orthogroups) > 1 else ("within_species_family" if len(species) == 1 else "cross_species_family")
        for query in query_ids:
            duplicate_context[query] = (f"DUP{duplicate_index:06d}", context, len(query_ids))
        duplicate_group_rows.append({
            "duplicate_group": f"DUP{duplicate_index:06d}", "sequence_md5": md5,
            "member_count": len(query_ids), "context": context,
            "orthogroups": ",".join(sorted(orthogroups)), "species": ",".join(sorted(species)),
            "query_ids": ",".join(sorted(query_ids)),
        })

    validation_rows = []
    risk_rows = []
    for row in proteins:
        query = row["query_id"]
        sequence = sequences[query]
        length = len(sequence)
        family = family_stats[row["orthogroup"]]
        length_ratio = length / family["median_length"] if family["median_length"] else 0.0
        hits = pfam_hits.get(query, [])
        direct_accessions = {hit["pfam_accession"] for hit in hits}
        direct_names = {hit["pfam_name"] for hit in hits}
        covered, coverage = interval_coverage([(hit["alignment_from"], hit["alignment_to"]) for hit in hits], length)
        architecture = ";".join(f"{hit['pfam_accession']}:{hit['alignment_from']}-{hit['alignment_to']}" for hit in hits)
        core = family["core_domains"]
        union = core | direct_accessions
        jaccard = len(core & direct_accessions) / len(union) if union else 1.0
        invalid = sorted(set(sequence) - VALID_AA)
        internal_stop = "*" in sequence[:-1]
        terminal_stop = sequence.endswith("*")
        unknown_fraction = sum(aa in {"X", "B", "Z", "J"} for aa in sequence) / length if length else 1.0
        entropy = sequence_entropy(sequence)
        taxid = parse_seed_taxid(row.get("seed_ortholog", ""))
        tax = taxa.get(taxid, {"taxid": taxid or "", "species": "", "rank": "", "track": "", "class": "unassigned"})
        risk_categories = flag_text(row, direct_names)
        flags = []
        severity = "pass"

        def add_flag(name, level):
            nonlocal severity
            flags.append(name)
            order = {"pass": 0, "review": 1, "review_high": 2, "exclude_recommended": 3}
            if order[level] > order[severity]:
                severity = level

        if length < 30:
            add_flag("very_short_lt30aa", "exclude_recommended")
        elif length < 50:
            add_flag("short_lt50aa", "review_high")
        elif length < 100:
            add_flag("short_lt100aa", "review")
        if invalid:
            add_flag("invalid_amino_acid_characters", "exclude_recommended")
        if internal_stop:
            add_flag("internal_stop_codon", "exclude_recommended")
        if terminal_stop:
            add_flag("terminal_stop_symbol", "review")
        if unknown_fraction > 0.05:
            add_flag("ambiguous_residue_fraction_gt5pct", "review")
        if entropy < 2.0 and length >= 50:
            add_flag("low_sequence_complexity", "review")
        if family["n"] >= 5 and length_ratio < 0.50:
            add_flag("severe_short_length_outlier", "review_high")
        elif family["n"] >= 5 and length_ratio < 0.67:
            add_flag("moderate_short_length_outlier", "review")
        if family["n"] >= 5 and length_ratio > 2.0:
            add_flag("severe_long_length_outlier", "review")
        duplicate = duplicate_context.get(query)
        if duplicate:
            if duplicate[1] == "cross_orthogroup":
                add_flag("exact_duplicate_cross_orthogroup", "review_high")
            elif duplicate[1] == "within_species_family":
                add_flag("exact_duplicate_within_species_family", "review")
            else:
                add_flag("exact_duplicate_cross_species_family", "review")
        if family["n_with_domains"] >= 3 and family["domain_rich_fraction"] >= 0.70 and not hits:
            add_flag("missing_domain_in_domain_rich_family", "review")
        if core and hits and jaccard < 0.25:
            add_flag("direct_domain_architecture_outlier", "review")
        if tax["class"] in {"bacteria", "archaea", "virus"}:
            add_flag(f"seed_taxonomy_{tax['class']}", "exclude_recommended")
        elif tax["class"] in {"fungi", "viridiplantae", "alveolata", "stramenopiles", "other_eukaryote"}:
            add_flag(f"non_metazoan_seed_{tax['class']}", "review_high")
        for category in risk_categories:
            add_flag(f"functional_risk_{category}", "review_high")

        transferred = {term for term in (row.get("PFAMs", "") or "").split(",") if term}
        domain_name_union = transferred | direct_names
        domain_name_intersection = transferred & direct_names
        transfer_jaccard = len(domain_name_intersection) / len(domain_name_union) if domain_name_union else 1.0
        dup_id, dup_kind, dup_n = duplicate if duplicate else ("", "", 1)
        out = {
            "group": args.group, "query_id": query, "species": row["species"], "gene_id": row["gene_id"],
            "orthogroup": row["orthogroup"], "validation_status": severity,
            "validation_flags": ";".join(sorted(set(flags))), "sequence_length_aa": length,
            "family_size": family["n"], "family_median_length_aa": f"{family['median_length']:.3f}",
            "length_ratio_to_family_median": f"{length_ratio:.6f}", "sequence_entropy_bits": f"{entropy:.6f}",
            "ambiguous_residue_fraction": f"{unknown_fraction:.6f}", "invalid_characters": "".join(invalid),
            "internal_stop": "yes" if internal_stop else "no", "terminal_stop": "yes" if terminal_stop else "no",
            "duplicate_group": dup_id, "duplicate_context": dup_kind, "duplicate_group_size": dup_n,
            "direct_pfam_domain_count": len(hits), "direct_pfam_unique_count": len(direct_accessions),
            "direct_pfam_covered_aa": covered, "direct_pfam_coverage_fraction": f"{coverage:.6f}",
            "direct_pfam_accessions": ",".join(sorted(direct_accessions)),
            "direct_pfam_names": ",".join(sorted(direct_names)), "direct_pfam_architecture": architecture,
            "family_core_pfam_accessions": ",".join(sorted(core)),
            "core_domain_jaccard": f"{jaccard:.6f}", "transferred_direct_pfam_jaccard": f"{transfer_jaccard:.6f}",
            "seed_ortholog": row.get("seed_ortholog", ""), "seed_taxid": tax.get("taxid", ""),
            "seed_species": tax.get("species", ""), "seed_lineage_class": tax.get("class", ""),
            "functional_risk_categories": ";".join(risk_categories),
            "description": row.get("Description", ""), "preferred_name": row.get("Preferred_name", ""),
            "transferred_pfams": row.get("PFAMs", ""), "go_terms": row.get("GOs", ""),
        }
        validation_rows.append(out)
        if risk_categories or tax["class"] not in {"metazoa", "unassigned"} or severity in {"review_high", "exclude_recommended"}:
            risk_rows.append(out)

    status_order = {"exclude_recommended": 0, "review_high": 1, "review": 2, "pass": 3}
    validation_rows.sort(key=lambda row: (status_order[row["validation_status"]], row["orthogroup"], row["species"], row["query_id"]))
    fields = list(validation_rows[0])
    write_tsv(args.output_dir / "sequence_validation.tsv", fields, validation_rows)
    write_tsv(args.output_dir / "high_risk_sequence_flags.tsv", fields, risk_rows)
    write_tsv(args.output_dir / "exact_duplicate_groups.tsv", [
        "duplicate_group", "sequence_md5", "member_count", "context", "orthogroups", "species", "query_ids"
    ], duplicate_group_rows)

    domain_fields = ["group", "species", "gene_id", "orthogroup"] + list(pfam_rows[0]) if pfam_rows else ["group", "query_id"]
    domain_out = []
    for hit in pfam_rows:
        meta = protein_by_query[hit["query_id"]]
        domain_out.append({"group": args.group, "species": meta["species"], "gene_id": meta["gene_id"], "orthogroup": meta["orthogroup"], **hit})
    write_tsv(args.output_dir / "direct_pfam_domains.tsv", domain_fields, domain_out)

    by_family_status = defaultdict(Counter)
    by_family_risk = defaultdict(Counter)
    by_family_tax = defaultdict(Counter)
    family_domain_outliers = Counter()
    for row in validation_rows:
        og = row["orthogroup"]
        by_family_status[og][row["validation_status"]] += 1
        for value in row["functional_risk_categories"].split(";"):
            if value:
                by_family_risk[og][value] += 1
        by_family_tax[og][row["seed_lineage_class"]] += 1
        if "direct_domain_architecture_outlier" in row["validation_flags"]:
            family_domain_outliers[og] += 1

    family_rows = []
    for orthogroup, stat in family_stats.items():
        meta = family_meta.get(orthogroup, {})
        statuses = by_family_status[orthogroup]
        risks = by_family_risk[orthogroup]
        priority = statuses["exclude_recommended"] * 100 + statuses["review_high"] * 10 + statuses["review"]
        family_rows.append({
            "group": args.group, "orthogroup": orthogroup, "family_size": stat["n"],
            "median_length_aa": f"{stat['median_length']:.3f}", "direct_pfam_sequences": stat["n_with_domains"],
            "direct_pfam_fraction": f"{stat['domain_rich_fraction']:.6f}",
            "family_core_pfam_accessions": ",".join(sorted(stat["core_domains"])),
            "pass_sequences": statuses["pass"], "review_sequences": statuses["review"],
            "review_high_sequences": statuses["review_high"], "exclude_recommended_sequences": statuses["exclude_recommended"],
            "domain_architecture_outliers": family_domain_outliers[orthogroup],
            "plastid_or_photosynthesis_flags": risks["plastid_or_photosynthesis"],
            "transposable_element_flags": risks["transposable_element"], "viral_or_phage_flags": risks["viral_or_phage"],
            "metazoan_seed_sequences": by_family_tax[orthogroup]["metazoa"],
            "non_metazoan_seed_sequences": sum(count for label, count in by_family_tax[orthogroup].items() if label not in {"metazoa", "unassigned"}),
            "validation_priority_score": priority, "step06_evidence_score": meta.get("evidence_score", ""),
            "best_event_node": meta.get("best_event_node", ""), "best_event_direction": meta.get("best_event_direction", ""),
            "best_event_change": meta.get("best_event_change", ""), "consensus_description": meta.get("consensus_description", ""),
        })
    family_rows.sort(key=lambda row: (-row["validation_priority_score"], row["orthogroup"]))
    write_tsv(args.output_dir / "family_validation_summary.tsv", list(family_rows[0]), family_rows)

    status_counts = Counter(row["validation_status"] for row in validation_rows)
    summary = [{
        "group": args.group, "total_sequences": len(validation_rows), "pass": status_counts["pass"],
        "review": status_counts["review"], "review_high": status_counts["review_high"],
        "exclude_recommended": status_counts["exclude_recommended"],
        "sequences_with_direct_pfam": sum(bool(row["direct_pfam_domain_count"]) for row in validation_rows),
        "exact_duplicate_groups": len(duplicate_group_rows),
        "non_metazoan_seed_flags": sum("non_metazoan_seed_" in row["validation_flags"] for row in validation_rows),
        "plastid_or_photosynthesis_flags": sum("plastid_or_photosynthesis" in row["functional_risk_categories"] for row in validation_rows),
        "transposable_element_flags": sum("transposable_element" in row["functional_risk_categories"] for row in validation_rows),
    }]
    write_tsv(args.output_dir / "validation_summary.tsv", list(summary[0]), summary)

    by_status = defaultdict(list)
    for row in validation_rows:
        by_status[row["validation_status"]].append(row["query_id"])
    write_fasta(args.output_dir / "validated_pass.faa", sequences, by_status["pass"])
    write_fasta(args.output_dir / "manual_review.faa", sequences, by_status["review"] + by_status["review_high"])
    write_fasta(args.output_dir / "exclude_recommended.faa", sequences, by_status["exclude_recommended"])
    provisional = [row["query_id"] for row in validation_rows if row["validation_status"] != "exclude_recommended"]
    write_fasta(args.output_dir / "provisional_retained.faa", sequences, provisional)
    shutil.copy2(args.domtblout, args.output_dir / "raw_hmmscan.domtblout")

    methods = f"""Step 9 Candidate Sequence Validation ({args.group})

Direct domain evidence
- HMMER hmmscan against the current Pfam-A profile library with model-specific curated GA thresholds (--cut_ga).
- Domain coordinates, coverage, architecture and family-level core domains are reported.

Sequence checks
- Protein length, invalid/ambiguous residues, internal stop symbols, sequence complexity and family-relative length outliers.
- Exact amino-acid duplicate detection within species/family, across species and across orthogroups.
- Direct Pfam architecture compared with family core domains and eggNOG-transferred Pfam labels.

Taxonomic and functional risk screening
- The taxonomic lineage of the eggNOG seed ortholog is classified using eggnog.taxa.db.
- Bacterial, archaeal or viral seed lineages are exclusion recommendations; non-metazoan eukaryotic seed lineages require high-priority review.
- Plastid/photosynthesis, transposable-element and viral/phage functional labels are flagged for review.

Triage policy
- pass: no triggered checks.
- review: moderate length, redundancy, low-complexity or domain-consistency concern.
- review_high: severe length/domain concern, non-metazoan eukaryotic seed, cross-orthogroup duplicate, or organelle/transposon/viral functional signal.
- exclude_recommended: hard sequence-format failure, internal stop, <30 aa, or bacterial/archaeal/viral seed lineage.
- Original inputs are never deleted. provisional_retained.faa excludes only exclude_recommended sequences; validated_pass.faa is the strict subset.

Limitations
- A non-metazoan seed hit is not proof of contamination and is never excluded on its own.
- Pfam architecture differences may reflect real domain gain/loss, alternative prediction or fragmentation; candidate families require manual review and, where possible, genome/transcript evidence.
- Taxonomic confirmation with a dedicated curated sequence database and genomic alignment remains advisable for the highest-priority candidates.
"""
    (args.output_dir / "METHODS_AND_LIMITATIONS.txt").write_text(methods)
    manifest = []
    for path in sorted(args.output_dir.iterdir()):
        if path.is_file() and path.name not in {"COMPLETE", "file_manifest.tsv"}:
            manifest.append({"file": path.name, "bytes": path.stat().st_size})
    write_tsv(args.output_dir / "file_manifest.tsv", ["file", "bytes"], manifest)
    (args.output_dir / "COMPLETE").write_text("complete\n")
    print(f"Validated {len(validation_rows)} {args.group} candidate proteins")


def cmd_combine(args):
    rows = []
    high_risk_families = []
    for group, directory in (("bivalvia", args.bivalvia), ("gastropoda", args.gastropoda)):
        if not (directory / "COMPLETE").exists():
            raise ValueError(f"Incomplete group validation: {directory}")
        rows.extend(read_tsv(directory / "validation_summary.tsv"))
        families = read_tsv(directory / "family_validation_summary.tsv")
        high_risk_families.extend([row for row in families if int(row["validation_priority_score"]) > 0])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "step09_summary.tsv", list(rows[0]), rows)
    high_risk_families.sort(key=lambda row: (-int(row["validation_priority_score"]), row["group"], row["orthogroup"]))
    if high_risk_families:
        write_tsv(args.output_dir / "cross_group_family_review_priorities.tsv", list(high_risk_families[0]), high_risk_families)
    else:
        write_tsv(args.output_dir / "cross_group_family_review_priorities.tsv", ["group", "orthogroup"], [])
    paths = [
        {"group": "bivalvia", "result_directory": str(args.bivalvia.resolve())},
        {"group": "gastropoda", "result_directory": str(args.gastropoda.resolve())},
    ]
    write_tsv(args.output_dir / "step09_result_directories.tsv", list(paths[0]), paths)
    (args.output_dir / "COMPLETE").write_text("complete\n")
    print("Combined Step 9 validation summaries")


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    validate = sub.add_parser("validate")
    validate.add_argument("--group", choices=["bivalvia", "gastropoda"], required=True)
    validate.add_argument("--fasta", type=Path, required=True)
    validate.add_argument("--annotations", type=Path, required=True)
    validate.add_argument("--family-table", type=Path, required=True)
    validate.add_argument("--domtblout", type=Path, required=True)
    validate.add_argument("--taxonomy-db", type=Path, required=True)
    validate.add_argument("--output-dir", type=Path, required=True)
    validate.set_defaults(func=cmd_validate)
    combine = sub.add_parser("combine")
    combine.add_argument("--bivalvia", type=Path, required=True)
    combine.add_argument("--gastropoda", type=Path, required=True)
    combine.add_argument("--output-dir", type=Path, required=True)
    combine.set_defaults(func=cmd_combine)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
