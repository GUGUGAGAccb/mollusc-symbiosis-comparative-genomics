#!/usr/bin/env python3
"""Summarise OrthoFinder orthogroup distributions without third-party packages."""

from __future__ import annotations

import argparse
import csv
import math
import os
import statistics
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


BIVALVIA = {
    "Acanthocardia_echinata",
    "Americardia_media",
    "Cerastoderma_edule",
    "Fragum_sueziense",
    "Tridacna_gigas",
    "Tridacna_maxima",
}

GASTROPODA = {
    "Aplysia_californica",
    "Berghia_stephanieae",
    "Bullacta_exarata",
    "Elysia_chlorotica",
    "Elysia_crispata",
    "Elysia_marginata",
    "Onchidella_celtica",
    "Plakobranchus_ocellatus",
}

GROUP_PATHS = {
    "combined": "orthofinder_results/core_complete80_v3_1_5/Results_Aug04",
    "bivalvia": "orthofinder_results/by_class/bivalvia_v3_1_5/Results_Aug04",
    "gastropoda": "orthofinder_results/by_class/gastropoda_v3_1_5/Results_Aug04",
}

EXPECTED_SPECIES = {
    "combined": BIVALVIA | GASTROPODA,
    "bivalvia": BIVALVIA,
    "gastropoda": GASTROPODA,
}

CATEGORY_ORDER = [
    "universal_single_copy",
    "universal_multicopy",
    "near_core",
    "shell",
    "species_specific",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_stats(path: Path) -> dict[str, str]:
    stats: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] and parts[0] not in stats:
                stats[parts[0]] = parts[1]
    return stats


def count_gene_ids(cell: str) -> int:
    if not cell.strip():
        return 0
    return sum(1 for item in cell.split(",") if item.strip())


def read_unassigned(path: Path, species: list[str]) -> dict[str, int]:
    counts = {name: 0 for name in species}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        for name in species:
            if name not in reader.fieldnames:
                raise ValueError(f"Species {name} missing from {path}")
        for row in reader:
            for name in species:
                counts[name] += count_gene_ids(row.get(name, ""))
    return counts


def classify(counts: list[int], near_core_fraction: float = 0.80) -> str:
    present = sum(value > 0 for value in counts)
    n_species = len(counts)
    near_core_minimum = math.ceil(n_species * near_core_fraction)
    if present == n_species and all(value == 1 for value in counts):
        return "universal_single_copy"
    if present == n_species:
        return "universal_multicopy"
    if present >= near_core_minimum:
        return "near_core"
    if present == 1:
        return "species_specific"
    return "shell"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def svg_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_bar_svg(path: Path, title: str, labels: list[str], values: list[float], suffix: str = "") -> None:
    width = 1100
    left = 300
    right = 90
    top = 90
    row_height = 48
    height = top + row_height * len(labels) + 70
    maximum = max(values) if values else 1.0
    plot_width = width - left - right
    colors = ["#2f75b5", "#4f9d69", "#d89b32", "#7d65a8", "#b85c5c"]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="40" y="45" font-family="Arial, sans-serif" font-size="25" font-weight="bold" fill="#17365d">{svg_escape(title)}</text>',
    ]
    for index, (label, value) in enumerate(zip(labels, values)):
        y = top + index * row_height
        bar_width = 0 if maximum == 0 else plot_width * value / maximum
        color = colors[index % len(colors)]
        lines.append(
            f'<text x="{left - 12}" y="{y + 25}" text-anchor="end" font-family="Arial, sans-serif" font-size="17" fill="#222222">{svg_escape(label)}</text>'
        )
        lines.append(
            f'<rect x="{left}" y="{y + 6}" width="{bar_width:.2f}" height="26" rx="3" fill="{color}"/>'
        )
        formatted = f"{value:,.1f}{suffix}" if suffix else f"{value:,.0f}"
        if bar_width > plot_width * 0.72:
            value_x = left + bar_width - 10
            anchor = "end"
            value_fill = "#ffffff"
        else:
            value_x = left + bar_width + 9
            anchor = "start"
            value_fill = "#222222"
        lines.append(
            f'<text x="{value_x:.2f}" y="{y + 25}" text-anchor="{anchor}" font-family="Arial, sans-serif" font-size="16" font-weight="bold" fill="{value_fill}">{formatted}</text>'
        )
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def analyse_group(project_root: Path, output_root: Path, group: str) -> dict[str, object]:
    base = project_root / GROUP_PATHS[group]
    gene_count_path = base / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    unassigned_path = base / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv"
    stats_path = base / "Comparative_Genomics_Statistics" / "Statistics_Overall.tsv"
    for path in (gene_count_path, unassigned_path, stats_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    group_output = output_root / group
    group_output.mkdir(parents=True, exist_ok=True)
    stats = read_stats(stats_path)
    rows: list[dict[str, object]] = []
    categories: Counter[str] = Counter()
    category_genes: Counter[str] = Counter()
    assigned: dict[str, int] = {}
    species: list[str] = []

    with gene_count_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) < 3 or header[0] != "Orthogroup" or header[-1] != "Total":
            raise ValueError(f"Unexpected GeneCount header in {gene_count_path}")
        species = header[1:-1]
        if set(species) != EXPECTED_SPECIES[group]:
            raise ValueError(f"Unexpected species for {group}: {species}")
        assigned = {name: 0 for name in species}
        for raw in reader:
            if not raw:
                continue
            orthogroup = raw[0]
            counts = [int(value) for value in raw[1:-1]]
            stated_total = int(raw[-1])
            calculated_total = sum(counts)
            if stated_total != calculated_total:
                raise ValueError(f"Total mismatch for {orthogroup}")
            for name, value in zip(species, counts):
                assigned[name] += value
            present_species = [name for name, value in zip(species, counts) if value > 0]
            category = classify(counts)
            categories[category] += 1
            category_genes[category] += calculated_total
            positive_counts = [value for value in counts if value > 0]
            row = {
                "orthogroup": orthogroup,
                "category": category,
                "n_species_present": len(present_species),
                "occupancy_fraction": f"{len(present_species) / len(species):.6f}",
                "total_copies": calculated_total,
                "mean_copies_all_species": f"{statistics.mean(counts):.6f}",
                "median_copies_all_species": f"{statistics.median(counts):.6f}",
                "max_copies": max(counts),
                "max_copy_species": ";".join(name for name, value in zip(species, counts) if value == max(counts)),
                "mean_copies_present_species": f"{statistics.mean(positive_counts):.6f}",
                "duplicated_in_any_species": "yes" if max(counts) > 1 else "no",
                "high_copy_ge100": "yes" if max(counts) >= 100 else "no",
                "present_species": ";".join(present_species),
            }
            for name, value in zip(species, counts):
                row[name] = value
            rows.append(row)

    distribution_fields = [
        "orthogroup",
        "category",
        "n_species_present",
        "occupancy_fraction",
        "total_copies",
        "mean_copies_all_species",
        "median_copies_all_species",
        "max_copies",
        "max_copy_species",
        "mean_copies_present_species",
        "duplicated_in_any_species",
        "high_copy_ge100",
        "present_species",
        *species,
    ]
    write_tsv(group_output / "orthogroup_distribution.tsv", distribution_fields, rows)

    total_orthogroups = len(rows)
    total_assigned = sum(assigned.values())
    summary_rows = []
    for category in CATEGORY_ORDER:
        count = categories[category]
        genes = category_genes[category]
        summary_rows.append(
            {
                "category": category,
                "orthogroups": count,
                "percent_orthogroups": f"{100 * count / total_orthogroups:.4f}",
                "assigned_genes": genes,
                "percent_assigned_genes": f"{100 * genes / total_assigned:.4f}",
            }
        )
    write_tsv(
        group_output / "category_summary.tsv",
        ["category", "orthogroups", "percent_orthogroups", "assigned_genes", "percent_assigned_genes"],
        summary_rows,
    )

    unassigned = read_unassigned(unassigned_path, species)
    species_rows = []
    for name in species:
        total = assigned[name] + unassigned[name]
        species_rows.append(
            {
                "species": name,
                "assigned_proteins": assigned[name],
                "unassigned_proteins": unassigned[name],
                "total_proteins": total,
                "assignment_rate_percent": f"{100 * assigned[name] / total:.4f}" if total else "0.0000",
            }
        )
    write_tsv(
        group_output / "species_assignment_summary.tsv",
        ["species", "assigned_proteins", "unassigned_proteins", "total_proteins", "assignment_rate_percent"],
        species_rows,
    )

    write_tsv(
        group_output / "species_specific_orthogroups.tsv",
        distribution_fields,
        [row for row in rows if row["category"] == "species_specific"],
    )
    write_tsv(
        group_output / "core_orthogroups.tsv",
        distribution_fields,
        [row for row in rows if row["category"] in {"universal_single_copy", "universal_multicopy", "near_core"}],
    )
    write_tsv(
        group_output / "high_copy_families_ge100.tsv",
        distribution_fields,
        [row for row in rows if row["high_copy_ge100"] == "yes"],
    )

    write_bar_svg(
        group_output / "orthogroup_categories.svg",
        f"{group}: orthogroup distribution",
        CATEGORY_ORDER,
        [float(categories[category]) for category in CATEGORY_ORDER],
    )
    write_bar_svg(
        group_output / "species_assignment_rates.svg",
        f"{group}: protein assignment rate",
        [row["species"] for row in species_rows],
        [float(row["assignment_rate_percent"]) for row in species_rows],
        "%",
    )

    expected_orthogroups = int(float(stats["Number of orthogroups"]))
    expected_assigned = int(float(stats["Number of genes in orthogroups"]))
    expected_genes = int(float(stats["Number of genes"]))
    observed_genes = total_assigned + sum(unassigned.values())
    checks = {
        "orthogroup_rows_match_statistics": total_orthogroups == expected_orthogroups,
        "assigned_gene_count_matches_statistics": total_assigned == expected_assigned,
        "total_gene_count_matches_statistics": observed_genes == expected_genes,
        "category_counts_sum_to_total": sum(categories.values()) == total_orthogroups,
    }
    if not all(checks.values()):
        raise RuntimeError(f"QC failed for {group}: {checks}")

    qc_rows = [{"check": key, "status": "PASS" if value else "FAIL"} for key, value in checks.items()]
    write_tsv(group_output / "qc_checks.tsv", ["check", "status"], qc_rows)

    return {
        "group": group,
        "base": base,
        "species": species,
        "rows": rows,
        "total_orthogroups": total_orthogroups,
        "total_genes": expected_genes,
        "assigned_genes": total_assigned,
        "assignment_rate": 100 * total_assigned / expected_genes,
        "categories": categories,
        "gene_count_path": gene_count_path,
        "unassigned_path": unassigned_path,
        "stats_path": stats_path,
    }


def write_combined_class_analysis(output_root: Path, combined: dict[str, object]) -> None:
    rows = combined["rows"]
    class_rows = []
    class_counts: Counter[str] = Counter()
    class_genes: Counter[str] = Counter()
    for row in rows:
        bivalvia_total = sum(int(row[name]) for name in BIVALVIA)
        gastropoda_total = sum(int(row[name]) for name in GASTROPODA)
        if bivalvia_total > 0 and gastropoda_total > 0:
            class_category = "shared_between_classes"
        elif bivalvia_total > 0:
            class_category = "bivalvia_only"
        elif gastropoda_total > 0:
            class_category = "gastropoda_only"
        else:
            raise RuntimeError(f"Empty orthogroup: {row['orthogroup']}")
        class_counts[class_category] += 1
        class_genes[class_category] += int(row["total_copies"])
        class_rows.append(
            {
                "orthogroup": row["orthogroup"],
                "class_category": class_category,
                "bivalvia_total_copies": bivalvia_total,
                "bivalvia_species_present": sum(int(row[name]) > 0 for name in BIVALVIA),
                "gastropoda_total_copies": gastropoda_total,
                "gastropoda_species_present": sum(int(row[name]) > 0 for name in GASTROPODA),
                "total_copies": row["total_copies"],
                "max_copies": row["max_copies"],
                "high_copy_ge100": row["high_copy_ge100"],
            }
        )
    fields = [
        "orthogroup",
        "class_category",
        "bivalvia_total_copies",
        "bivalvia_species_present",
        "gastropoda_total_copies",
        "gastropoda_species_present",
        "total_copies",
        "max_copies",
        "high_copy_ge100",
    ]
    write_tsv(output_root / "combined" / "class_distribution.tsv", fields, class_rows)
    categories = ["shared_between_classes", "bivalvia_only", "gastropoda_only"]
    summary = [
        {
            "class_category": category,
            "orthogroups": class_counts[category],
            "percent_orthogroups": f"{100 * class_counts[category] / len(class_rows):.4f}",
            "assigned_genes": class_genes[category],
        }
        for category in categories
    ]
    write_tsv(
        output_root / "combined" / "class_distribution_summary.tsv",
        ["class_category", "orthogroups", "percent_orthogroups", "assigned_genes"],
        summary,
    )
    for category in ("bivalvia_only", "gastropoda_only"):
        write_tsv(
            output_root / "combined" / f"{category}_orthogroups.tsv",
            fields,
            [row for row in class_rows if row["class_category"] == category],
        )
    write_bar_svg(
        output_root / "combined" / "class_distribution.svg",
        "Combined analysis: orthogroups by class distribution",
        categories,
        [float(class_counts[category]) for category in categories],
    )


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    output_root = args.output.resolve()
    output_root.mkdir(parents=True, exist_ok=False)

    analyses = {group: analyse_group(project_root, output_root, group) for group in GROUP_PATHS}
    write_combined_class_analysis(output_root, analyses["combined"])

    overall_rows = []
    for group in ("combined", "bivalvia", "gastropoda"):
        result = analyses[group]
        categories = result["categories"]
        overall_rows.append(
            {
                "group": group,
                "species": len(result["species"]),
                "proteins": result["total_genes"],
                "assigned_proteins": result["assigned_genes"],
                "assignment_rate_percent": f"{result['assignment_rate']:.4f}",
                "orthogroups": result["total_orthogroups"],
                "universal_single_copy": categories["universal_single_copy"],
                "universal_multicopy": categories["universal_multicopy"],
                "near_core": categories["near_core"],
                "shell": categories["shell"],
                "species_specific": categories["species_specific"],
            }
        )
    write_tsv(
        output_root / "step04_summary.tsv",
        [
            "group",
            "species",
            "proteins",
            "assigned_proteins",
            "assignment_rate_percent",
            "orthogroups",
            "universal_single_copy",
            "universal_multicopy",
            "near_core",
            "shell",
            "species_specific",
        ],
        overall_rows,
    )

    manifest_rows = []
    for group, result in analyses.items():
        for role in ("gene_count_path", "unassigned_path", "stats_path"):
            path = result[role]
            stat = path.stat()
            manifest_rows.append(
                {
                    "group": group,
                    "role": role,
                    "path": str(path),
                    "size_bytes": stat.st_size,
                    "mtime_utc": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
                }
            )
    write_tsv(
        output_root / "input_manifest.tsv",
        ["group", "role", "path", "size_bytes", "mtime_utc"],
        manifest_rows,
    )

    readme = [
        "# Step 04 — Orthogroup distribution analysis",
        "",
        "Status: complete",
        "",
        "Definitions:",
        "",
        "- universal_single_copy: exactly one copy in every species",
        "- universal_multicopy: present in every species and duplicated in at least one species",
        "- near_core: present in at least 80% but not all species",
        "- shell: present in two or more species but below the near-core threshold",
        "- species_specific: present in exactly one species",
        "- bivalvia_only / gastropoda_only: present only in that class in the combined 14-species OrthoFinder run",
        "",
        "Quality control:",
        "",
        "- GeneCount row totals were validated.",
        "- Orthogroup, assigned-gene and total-gene counts were matched to OrthoFinder Statistics_Overall.tsv.",
        "- Unassigned protein counts were calculated from Orthogroups_UnassignedGenes.tsv.",
        "",
        "Primary outputs:",
        "",
        "- step04_summary.tsv",
        "- <group>/orthogroup_distribution.tsv",
        "- <group>/category_summary.tsv",
        "- <group>/species_assignment_summary.tsv",
        "- combined/class_distribution.tsv",
        "- combined/{bivalvia_only,gastropoda_only}_orthogroups.tsv",
        "- SVG summary figures and qc_checks.tsv",
        "",
        f"Generated: {datetime.now(tz=timezone.utc).isoformat()}",
    ]
    (output_root / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    (output_root / "COMPLETE").write_text(datetime.now(tz=timezone.utc).isoformat() + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
