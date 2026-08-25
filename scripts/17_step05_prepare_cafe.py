#!/usr/bin/env python3
"""Root an IQ-TREE topology, make a relative ultrametric tree, and prepare CAFE5 input."""

from __future__ import annotations

import argparse
import csv
import math
from datetime import datetime, timezone
from pathlib import Path

from ete4 import Tree


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--group", required=True, choices=("bivalvia", "gastropoda"))
    parser.add_argument("--gene-count", required=True, type=Path)
    parser.add_argument("--iqtree-tree", required=True, type=Path)
    parser.add_argument("--reference-rooted-tree", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--root-age", type=float, default=100.0)
    parser.add_argument("--max-family-size", type=int, default=99)
    return parser.parse_args()


def leaves(node) -> set[str]:
    # ETE4 removed ETE3's get_leaf_names() helper; leaves() is stable in ETE4.
    return {leaf.name for leaf in node.leaves()}


def root_to_node_distance(node) -> float:
    distance = 0.0
    current = node
    while current.up is not None:
        distance += float(current.dist)
        current = current.up
    return distance


def root_like_reference(target: Tree, reference: Tree) -> tuple[set[str], set[str]]:
    target_leaves = leaves(target)
    reference_leaves = leaves(reference)
    if target_leaves != reference_leaves:
        missing = sorted(reference_leaves - target_leaves)
        extra = sorted(target_leaves - reference_leaves)
        raise ValueError(f"Tree leaf mismatch; missing={missing}; extra={extra}")
    reference_children = list(reference.children)
    if len(reference_children) != 2:
        raise ValueError("Reference species tree is not rooted and binary at the root")
    reference_sides = [leaves(child) for child in reference_children]
    if not all(reference_sides):
        raise ValueError("Empty root partition in reference tree")

    current_children = list(target.children)
    if len(current_children) == 2:
        current_sides = [leaves(child) for child in current_children]
        if {frozenset(side) for side in current_sides} == {frozenset(side) for side in reference_sides}:
            return current_sides[0], current_sides[1]

    preferred = min(reference_sides, key=lambda side: (len(side), sorted(side)))
    alternatives = [preferred, target_leaves - preferred]
    outgroup_node = None
    for wanted in alternatives:
        for node in target.traverse():
            if node is target:
                continue
            if leaves(node) == wanted:
                outgroup_node = node
                break
        if outgroup_node is not None:
            break
    if outgroup_node is None:
        raise ValueError("IQ-TREE topology does not contain the OrthoFinder root bipartition")

    target.set_outgroup(outgroup_node)
    rooted_children = list(target.children)
    if len(rooted_children) != 2:
        raise ValueError("Rerooted tree is not binary at the root")
    rooted_sides = [leaves(child) for child in rooted_children]
    if {frozenset(side) for side in rooted_sides} != {frozenset(side) for side in reference_sides}:
        raise ValueError("Rerooted IQ-TREE does not match the reference root bipartition")
    return rooted_sides[0], rooted_sides[1]


def make_relative_ultrametric(tree: Tree, root_age: float) -> tuple[float, float]:
    if root_age <= 0:
        raise ValueError("root_age must be positive")
    leaf_nodes = list(tree.leaves())
    heights = [root_to_node_distance(leaf) for leaf in leaf_nodes]
    maximum = max(heights)
    if maximum <= 0:
        raise ValueError("Tree has no positive branch lengths")

    # First extend terminal branches only enough to reach the longest observed root-to-tip path.
    for leaf, height in zip(leaf_nodes, heights):
        leaf.dist = float(leaf.dist) + (maximum - height)

    # Uniform scaling changes only the time unit and preserves relative internal branch lengths.
    scale = root_age / maximum
    for node in tree.traverse():
        if node.up is not None:
            node.dist = float(node.dist) * scale

    final_heights = [root_to_node_distance(leaf) for leaf in tree.leaves()]
    spread = max(final_heights) - min(final_heights)
    if spread > max(1e-7, root_age * 1e-8):
        raise RuntimeError(f"Ultrametric validation failed; root-to-tip spread={spread}")
    if any(float(node.dist) < 0 for node in tree.traverse() if node.up is not None):
        raise RuntimeError("Negative branch length after ultrametric conversion")
    return maximum, scale


def prepare_families(
    gene_count_path: Path,
    species_order: list[str],
    root_side_a: set[str],
    root_side_b: set[str],
    max_family_size: int,
    output: Path,
) -> dict[str, int]:
    primary_path = output / "cafe_input_primary.tsv"
    large_path = output / "cafe_input_large_ge100.tsv"
    root_absent_path = output / "cafe_input_root_absent.tsv"
    header = ["Desc", "Family ID", *species_order]
    counts = {"total": 0, "primary": 0, "large": 0, "root_absent": 0}

    with gene_count_path.open(encoding="utf-8", newline="") as source, \
        primary_path.open("w", encoding="utf-8", newline="") as primary_handle, \
        large_path.open("w", encoding="utf-8", newline="") as large_handle, \
        root_absent_path.open("w", encoding="utf-8", newline="") as absent_handle:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing GeneCount header: {gene_count_path}")
        source_species = [name for name in reader.fieldnames if name not in {"Orthogroup", "Total"}]
        if set(source_species) != set(species_order):
            raise ValueError("GeneCount species do not match tree leaves")
        writers = {
            "primary": csv.writer(primary_handle, delimiter="\t", lineterminator="\n"),
            "large": csv.writer(large_handle, delimiter="\t", lineterminator="\n"),
            "root_absent": csv.writer(absent_handle, delimiter="\t", lineterminator="\n"),
        }
        for writer in writers.values():
            writer.writerow(header)
        for row in reader:
            counts["total"] += 1
            family = row["Orthogroup"]
            values = {name: int(row[name]) for name in species_order}
            present_a = any(values[name] > 0 for name in root_side_a)
            present_b = any(values[name] > 0 for name in root_side_b)
            output_row = ["(null)", family, *(values[name] for name in species_order)]
            if not (present_a and present_b):
                category = "root_absent"
            elif max(values.values()) > max_family_size:
                category = "large"
            else:
                category = "primary"
            writers[category].writerow(output_row)
            counts[category] += 1
    if counts["primary"] == 0:
        raise RuntimeError("No root-present gene families remain after filtering")
    if counts["total"] != counts["primary"] + counts["large"] + counts["root_absent"]:
        raise RuntimeError("Family partition counts do not sum to the input total")
    return counts


def main() -> None:
    args = parse_args()
    for path in (args.gene_count, args.iqtree_tree, args.reference_rooted_tree):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)

    reference = Tree(str(args.reference_rooted_tree), parser=1)
    target = Tree(str(args.iqtree_tree), parser=1)
    root_side_a, root_side_b = root_like_reference(target, reference)
    original_max_height, scale = make_relative_ultrametric(target, args.root_age)

    rooted_tree_path = output / "SpeciesTree_IQTREE.rooted_relative_ultrametric.tre"
    target.write(outfile=str(rooted_tree_path), parser=5)
    species_order = sorted(leaves(target))
    family_counts = prepare_families(
        args.gene_count,
        species_order,
        root_side_a,
        root_side_b,
        args.max_family_size,
        output,
    )

    with (output / "cafe_input_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["group", args.group])
        writer.writerow(["species", len(species_order)])
        writer.writerow(["root_side_a", ";".join(sorted(root_side_a))])
        writer.writerow(["root_side_b", ";".join(sorted(root_side_b))])
        writer.writerow(["relative_root_age", args.root_age])
        writer.writerow(["original_max_root_to_tip", f"{original_max_height:.12g}"])
        writer.writerow(["uniform_scale_factor", f"{scale:.12g}"])
        writer.writerow(["max_family_size_primary", args.max_family_size])
        for key in ("total", "primary", "large", "root_absent"):
            writer.writerow([f"families_{key}", family_counts[key]])

    methods = [
        f"Group: {args.group}",
        "IQ-TREE topology was rooted using the root bipartition in OrthoFinder SpeciesTree_rooted.txt.",
        "Terminal branches were extended to the longest observed root-to-tip distance, then all branches were uniformly scaled to a relative root age of 100.",
        "The tree is ultrametric in relative units; it is not fossil- or TimeTree-calibrated.",
        "Families absent from either side of the root split were excluded from the primary CAFE5 likelihood analysis.",
        f"Families with any species count greater than {args.max_family_size} were placed in a separate large-family table.",
        "Primary CAFE5 results are exploratory screens until calibration and alternative ultrametric methods are assessed.",
        f"Generated: {datetime.now(tz=timezone.utc).isoformat()}",
    ]
    (output / "METHODS_AND_LIMITATIONS.txt").write_text("\n".join(methods) + "\n", encoding="utf-8")
    (output / "PREPARE_COMPLETE").write_text(datetime.now(tz=timezone.utc).isoformat() + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
