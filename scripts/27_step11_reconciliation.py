#!/usr/bin/env python3
"""Step 11: duplication-loss LCA reconciliation of candidate gene trees."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade


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
                    raise ValueError(f"Sequence before FASTA header: {path}")
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


def fmt(value):
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def parse_support(name, confidence=None):
    sh_alrt = None
    ufboot = None
    if name:
        match = re.fullmatch(r"\s*([0-9.]+)\s*/\s*([0-9.]+)\s*", str(name))
        if match:
            sh_alrt = number(match.group(1), None)
            ufboot = number(match.group(2), None)
        else:
            single = re.fullmatch(r"\s*([0-9.]+)\s*", str(name))
            if single:
                ufboot = number(single.group(1), None)
    if ufboot is None and confidence is not None:
        ufboot = number(confidence, None)
    if sh_alrt is not None and ufboot is not None and sh_alrt >= 80 and ufboot >= 95:
        support_class = "strong"
    elif (sh_alrt is not None and sh_alrt >= 80) or (ufboot is not None and ufboot >= 95):
        support_class = "partial_strong"
    elif (sh_alrt is not None and sh_alrt >= 60) or (ufboot is not None and ufboot >= 70):
        support_class = "moderate"
    elif sh_alrt is None and ufboot is None:
        support_class = "unknown"
    else:
        support_class = "low"
    return sh_alrt, ufboot, support_class


class SpeciesModel:
    def __init__(self, tree_path: Path):
        self.tree_path = tree_path
        self.tree = Phylo.read(str(tree_path), "newick")
        self.root = self.tree.root
        self.parent = {self.root: None}
        self.depth = {self.root: 0}
        self.children = {}
        self.descendants = {}
        self.node_id = {}
        self.by_species = {}
        self._walk(self.root)
        for terminal in self.tree.get_terminals():
            if not terminal.name:
                raise ValueError(f"Unnamed species-tree terminal in {tree_path}")
            if terminal.name in self.by_species:
                raise ValueError(f"Duplicate species-tree terminal: {terminal.name}")
            self.by_species[terminal.name] = terminal
        internal = []
        for node in self.tree.find_clades(order="postorder"):
            leaves = frozenset(term.name for term in node.get_terminals())
            self.descendants[node] = leaves
            if not node.is_terminal():
                internal.append(node)
        for idx, node in enumerate(sorted(internal, key=lambda n: (len(self.descendants[n]), sorted(self.descendants[n]))), 1):
            self.node_id[node] = "SROOT" if node is self.root else f"S{idx:03d}"
        for terminal in self.tree.get_terminals():
            self.node_id[terminal] = terminal.name

    def _walk(self, node):
        self.children[node] = list(node.clades)
        for child in node.clades:
            self.parent[child] = node
            self.depth[child] = self.depth[node] + 1
            self._walk(child)

    def lca_nodes(self, nodes):
        nodes = list(nodes)
        if not nodes:
            raise ValueError("Cannot calculate LCA of an empty node set")
        current = nodes[0]
        for other in nodes[1:]:
            a, b = current, other
            while self.depth[a] > self.depth[b]:
                a = self.parent[a]
            while self.depth[b] > self.depth[a]:
                b = self.parent[b]
            while a is not b:
                a = self.parent[a]
                b = self.parent[b]
            current = a
        return current

    def lca_species(self, species):
        missing = sorted(set(species) - set(self.by_species))
        if missing:
            raise ValueError(f"Gene tree contains species absent from species tree: {missing}")
        return self.lca_nodes(self.by_species[name] for name in species)

    def distance(self, ancestor, descendant):
        if self.depth[descendant] < self.depth[ancestor]:
            raise ValueError("Species mapping is not ancestor-descendant")
        distance = 0
        current = descendant
        while current is not ancestor:
            current = self.parent[current]
            distance += 1
            if current is None:
                raise ValueError("Species mapping is not ancestor-descendant")
        return distance

    def path(self, ancestor, descendant):
        reverse = []
        current = descendant
        while current is not ancestor:
            reverse.append(current)
            current = self.parent[current]
            if current is None:
                raise ValueError("Species mapping is not ancestor-descendant")
        return list(reversed(reverse))

    def clade_text(self, node):
        return ";".join(sorted(self.descendants[node]))

    def branch_type(self, node):
        if node is self.root:
            return "root"
        return "terminal" if node.is_terminal() else "internal"

    def branch_rows(self, group):
        rows = []
        for node in self.tree.find_clades(order="preorder"):
            rows.append({
                "group": group,
                "species_node": self.node_id[node],
                "parent_species_node": self.node_id[self.parent[node]] if self.parent[node] is not None else "",
                "branch_type": self.branch_type(node),
                "descendant_species_count": len(self.descendants[node]),
                "descendant_species": self.clade_text(node),
            })
        return rows


@dataclass
class OrientedNode:
    original: object | None
    original_id: int | None
    incoming_length: float | None
    children: list["OrientedNode"] = field(default_factory=list)
    leaf_name: str = ""
    species: frozenset[str] = frozenset()
    species_node: object | None = None
    event: str = ""
    losses: int = 0
    gene_node_id: str = ""

    @property
    def is_leaf(self):
        return not self.children


def build_gene_graph(tree):
    nodes = list(tree.find_clades(order="preorder"))
    ids = {node: idx for idx, node in enumerate(nodes)}
    adjacency = {node: [] for node in nodes}
    for parent in nodes:
        for child in parent.clades:
            length = child.branch_length if child.branch_length is not None else 0.0
            adjacency[parent].append((child, length))
            adjacency[child].append((parent, length))
    edges = []
    for left in nodes:
        for right, length in adjacency[left]:
            if ids[left] < ids[right]:
                edges.append((left, right, length))
    return nodes, ids, adjacency, edges


def build_side(node, parent, incoming_length, adjacency, ids):
    children = [
        build_side(neighbour, node, length, adjacency, ids)
        for neighbour, length in adjacency[node]
        if neighbour is not parent
    ]
    leaf_name = node.name if not children else ""
    return OrientedNode(node, ids[node], incoming_length, children, leaf_name or "")


def orient_at_edge(edge, adjacency, ids):
    left, right, length = edge
    half = (length or 0.0) / 2.0
    return OrientedNode(None, None, None, [
        build_side(left, right, half, adjacency, ids),
        build_side(right, left, half, adjacency, ids),
    ])


def postorder(node):
    for child in node.children:
        yield from postorder(child)
    yield node


def preorder(node):
    yield node
    for child in node.children:
        yield from preorder(child)


def descendant_leaf_names(node):
    if node.is_leaf:
        return frozenset([node.leaf_name])
    leaves = set()
    for child in node.children:
        leaves.update(descendant_leaf_names(child))
    return frozenset(leaves)


def root_signature(root):
    sides = [tuple(sorted(descendant_leaf_names(child))) for child in root.children]
    sides.sort(key=lambda side: (len(side), side))
    return "|".join(sides[0]) + "||" + "|".join(sides[1])


def evaluate_reconciliation(root, label_to_species, species_model):
    duplications = 0
    speciations = 0
    losses = 0
    for node in postorder(root):
        if node.is_leaf:
            if node.leaf_name not in label_to_species:
                raise ValueError(f"Gene-tree label absent from sequence map: {node.leaf_name}")
            species = label_to_species[node.leaf_name]
            node.species = frozenset([species])
            node.species_node = species_model.by_species[species]
            node.event = "leaf"
            continue
        merged = set()
        for child in node.children:
            merged.update(child.species)
        node.species = frozenset(merged)
        node.species_node = species_model.lca_species(merged)
        is_duplication = any(child.species_node is node.species_node for child in node.children)
        node.event = "duplication" if is_duplication else "speciation"
        if is_duplication:
            duplications += 1
        else:
            speciations += 1
        node.losses = 0
        for child in node.children:
            distance = species_model.distance(node.species_node, child.species_node)
            edge_losses = max(0, distance - (0 if is_duplication else 1))
            node.losses += edge_losses
            losses += edge_losses
    return duplications, losses, speciations


def assign_gene_node_ids(root):
    counter = 0
    for node in postorder(root):
        if not node.is_leaf:
            counter += 1
            node.gene_node_id = f"G{counter:04d}"


def support_for_node(node):
    if node.original is None:
        return None, None, "root_induced"
    return parse_support(node.original.name, node.original.confidence)


def loss_events_for_edge(parent, child, species_model):
    path = species_model.path(parent.species_node, child.species_node)
    start = 0 if parent.event == "duplication" else 1
    events = []
    current = parent.species_node
    for idx, followed in enumerate(path):
        if idx >= start:
            siblings = [candidate for candidate in species_model.children[current] if candidate is not followed]
            for sibling in siblings:
                events.append(sibling)
        current = followed
    return events


def safe_newick_label(text):
    return re.sub(r"[^A-Za-z0-9_.@=-]+", "_", text)


def restore_identical_leaves(gene_tree, expected_labels, fasta_path):
    """Reinsert IQ-TREE-pruned identical sequences as zero-length polytomies."""
    present = {terminal.name for terminal in gene_tree.get_terminals()}
    missing = sorted(set(expected_labels) - present)
    if not missing:
        return 0
    sequences = read_fasta(fasta_path)
    absent_fasta = sorted((set(expected_labels) | present) - set(sequences))
    if absent_fasta:
        raise ValueError(f"Tree/map labels absent from input FASTA: {absent_fasta[:10]}")
    sequence_to_present = defaultdict(list)
    for label in present:
        sequence_to_present[sequences[label]].append(label)
    anchor_to_missing = defaultdict(list)
    unresolved = []
    for label in missing:
        anchors = sorted(sequence_to_present.get(sequences[label], []))
        if not anchors:
            unresolved.append(label)
        else:
            anchor_to_missing[anchors[0]].append(label)
    if unresolved:
        raise ValueError(
            f"IQ-TREE omitted {len(unresolved)} non-identical/unrecoverable labels: {unresolved[:10]}"
        )
    parent = {gene_tree.root: None}
    terminal_by_name = {}
    for node in gene_tree.find_clades(order="preorder"):
        if node.is_terminal():
            terminal_by_name[node.name] = node
        for child in node.clades:
            parent[child] = node
    for anchor_name, labels in sorted(anchor_to_missing.items()):
        anchor = terminal_by_name[anchor_name]
        old_parent = parent[anchor]
        old_length = anchor.branch_length
        anchor.branch_length = 0.0
        restored = [Clade(branch_length=0.0, name=label) for label in sorted(labels)]
        polytomy = Clade(branch_length=old_length, clades=[anchor, *restored])
        if old_parent is None:
            gene_tree.root = polytomy
        else:
            index = old_parent.clades.index(anchor)
            old_parent.clades[index] = polytomy
    return len(missing)


def serialize_reconciled(node, species_model):
    if node.is_leaf:
        body = safe_newick_label(node.leaf_name)
    else:
        children = ",".join(serialize_reconciled(child, species_model) for child in node.children)
        event = "D" if node.event == "duplication" else "S"
        sh_alrt, ufboot, _ = support_for_node(node)
        support = ""
        if sh_alrt is not None or ufboot is not None:
            support = f"_a{fmt(sh_alrt) or 'NA'}_b{fmt(ufboot) or 'NA'}"
        label = f"{event}_{node.gene_node_id}@{species_model.node_id[node.species_node]}{support}"
        body = f"({children}){safe_newick_label(label)}"
    if node.incoming_length is not None:
        body += f":{node.incoming_length:.10g}"
    return body


def reconcile_family(group, selected, tree_row, species_model, outdir, duplication_cost, loss_cost):
    orthogroup = selected["orthogroup"]
    map_rows = read_tsv(Path(selected["sequence_map"]))
    full_label_to_species = {row["tree_label"]: row["species"] for row in map_rows}
    map_by_label = {row["tree_label"]: row for row in map_rows}
    tree_path = Path(tree_row["treefile"])
    if not tree_path.exists():
        raise FileNotFoundError(f"Missing gene tree: {tree_path}")
    gene_tree = Phylo.read(str(tree_path), "newick")
    reconciliation_alignment = tree_path.parent / "trimmed.faa"
    if not reconciliation_alignment.exists():
        reconciliation_alignment = Path(selected["input_fasta"])
    alignment_labels = set(read_fasta(reconciliation_alignment))
    unexpected_alignment = sorted(alignment_labels - set(full_label_to_species))
    if unexpected_alignment:
        raise ValueError(f"Alignment contains labels absent from sequence map: {unexpected_alignment[:10]}")
    alignment_filtered = sorted(set(full_label_to_species) - alignment_labels)
    label_to_species = {
        label: species for label, species in full_label_to_species.items() if label in alignment_labels
    }
    reinserted_identical = restore_identical_leaves(
        gene_tree, label_to_species, reconciliation_alignment
    )
    terminal_names = [terminal.name for terminal in gene_tree.get_terminals()]
    if len(terminal_names) != len(set(terminal_names)):
        raise ValueError(f"Duplicate terminal labels in {orthogroup}")
    if set(terminal_names) != set(label_to_species):
        raise ValueError(
            f"Tree/map mismatch for {orthogroup}: tree_only={len(set(terminal_names)-set(label_to_species))}, "
            f"map_only={len(set(label_to_species)-set(terminal_names))}"
        )
    nodes, ids, adjacency, edges = build_gene_graph(gene_tree)
    candidates = []
    for edge in edges:
        root = orient_at_edge(edge, adjacency, ids)
        duplications, losses, speciations = evaluate_reconciliation(root, label_to_species, species_model)
        signature = root_signature(root)
        score = duplication_cost * duplications + loss_cost * losses
        candidates.append({
            "root": root,
            "duplications": duplications,
            "losses": losses,
            "speciations": speciations,
            "score": score,
            "signature": signature,
        })
    candidates.sort(key=lambda item: (item["score"], item["duplications"], item["losses"], item["signature"]))
    chosen = candidates[0]
    optimal = [item for item in candidates if abs(item["score"] - chosen["score"]) < 1e-9]
    root = chosen["root"]
    assign_gene_node_ids(root)

    duplication_frequency = Counter()
    for candidate in optimal:
        for node in preorder(candidate["root"]):
            if node.event != "duplication":
                continue
            key = (node.original_id, species_model.node_id[node.species_node])
            duplication_frequency[key] += 1

    duplication_rows = []
    loss_rows = []
    loss_counter = 0
    strong_duplications = 0
    low_support_duplications = 0
    for node in preorder(root):
        if node.event == "duplication":
            sh_alrt, ufboot, support_class = support_for_node(node)
            robustness = duplication_frequency[(node.original_id, species_model.node_id[node.species_node])] / len(optimal)
            if support_class == "strong":
                strong_duplications += 1
            if support_class in {"low", "unknown", "root_induced"}:
                low_support_duplications += 1
            leaves = sorted(descendant_leaf_names(node))
            duplication_rows.append({
                "group": group,
                "orthogroup": orthogroup,
                "gene_node_id": node.gene_node_id,
                "species_node": species_model.node_id[node.species_node],
                "species_branch_type": species_model.branch_type(node.species_node),
                "mapped_clade_species_count": len(species_model.descendants[node.species_node]),
                "mapped_clade_species": species_model.clade_text(node.species_node),
                "descendant_gene_count": len(leaves),
                "descendant_species_count": len(node.species),
                "descendant_species": ";".join(sorted(node.species)),
                "descendant_tree_labels": ";".join(leaves),
                "sh_alrt": fmt(sh_alrt),
                "ufboot": fmt(ufboot),
                "support_class": support_class,
                "optimal_root_robustness": f"{robustness:.6f}",
                "root_induced_event": "yes" if node.original is None else "no",
                "cafe_best_event_node": selected.get("best_event_node", ""),
                "cafe_best_event_direction": selected.get("best_event_direction", ""),
                "cafe_best_event_change": selected.get("best_event_change", ""),
                "selection_score": selected.get("selection_score", ""),
            })
        if node.is_leaf:
            continue
        for child in node.children:
            for lost_branch in loss_events_for_edge(node, child, species_model):
                loss_counter += 1
                sh_alrt, ufboot, support_class = support_for_node(child)
                loss_rows.append({
                    "group": group,
                    "orthogroup": orthogroup,
                    "loss_event_id": f"L{loss_counter:04d}",
                    "parent_gene_node_id": node.gene_node_id,
                    "child_gene_node_id": child.gene_node_id if not child.is_leaf else child.leaf_name,
                    "lost_species_node": species_model.node_id[lost_branch],
                    "lost_branch_type": species_model.branch_type(lost_branch),
                    "lost_clade_species_count": len(species_model.descendants[lost_branch]),
                    "lost_clade_species": species_model.clade_text(lost_branch),
                    "edge_sh_alrt": fmt(sh_alrt),
                    "edge_ufboot": fmt(ufboot),
                    "edge_support_class": support_class,
                })

    family_dir = outdir / "reconciled_trees" / orthogroup
    family_dir.mkdir(parents=True, exist_ok=True)
    reconciled_tree_path = family_dir / f"{orthogroup}.reconciled.nwk"
    reconciled_tree_path.write_text(serialize_reconciled(root, species_model) + ";\n")
    root_candidates = family_dir / f"{orthogroup}.root_candidates.tsv"
    candidate_rows = [{
        "orthogroup": orthogroup,
        "rank": rank,
        "root_signature": item["signature"],
        "duplications": item["duplications"],
        "losses": item["losses"],
        "speciations": item["speciations"],
        "weighted_dl_cost": f"{item['score']:.6f}",
        "optimal": "yes" if item in optimal else "no",
    } for rank, item in enumerate(candidates, 1)]
    write_tsv(root_candidates, list(candidate_rows[0]), candidate_rows)

    family_summary = {
        "group": group,
        "rank": selected.get("rank", ""),
        "orthogroup": orthogroup,
        "input_sequences": selected.get("n_sequences", ""),
        "input_species": selected.get("n_species", ""),
        "alignment_filtered_sequences": len(alignment_filtered),
        "reinserted_iqtree_identical_sequences": reinserted_identical,
        "reconciled_sequences": len(label_to_species),
        "reconciled_species": len(set(label_to_species.values())),
        "best_model": tree_row.get("best_model", ""),
        "cafe_best_event_node": selected.get("best_event_node", ""),
        "cafe_best_event_direction": selected.get("best_event_direction", ""),
        "cafe_best_event_change": selected.get("best_event_change", ""),
        "selection_score": selected.get("selection_score", ""),
        "chosen_root_signature": chosen["signature"],
        "tested_root_edges": len(candidates),
        "equally_optimal_roots": len(optimal),
        "duplications": chosen["duplications"],
        "strong_support_duplications": strong_duplications,
        "low_or_unknown_support_duplications": low_support_duplications,
        "losses": chosen["losses"],
        "speciations": chosen["speciations"],
        "weighted_dl_cost": f"{chosen['score']:.6f}",
        "duplication_species_nodes": ";".join(sorted({row["species_node"] for row in duplication_rows})),
        "loss_species_nodes": ";".join(sorted({row["lost_species_node"] for row in loss_rows})),
        "reconciled_tree": str(reconciled_tree_path),
        "root_candidates_table": str(root_candidates),
    }
    filtered_rows = [{
        "group": group,
        "orthogroup": orthogroup,
        "tree_label": label,
        "species": map_by_label[label].get("species", ""),
        "gene_id": map_by_label[label].get("gene_id", ""),
        "reason": "removed_by_alignment_trimming_before_iqtree",
    } for label in alignment_filtered]
    return family_summary, duplication_rows, loss_rows, filtered_rows


def cmd_reconcile_group(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    step10 = args.step10.resolve()
    selected_rows = read_tsv(step10 / f"selected_{args.group}.tsv")
    tree_rows = {row["orthogroup"]: row for row in read_tsv(step10 / "gene_tree_summary.tsv") if row["group"] == args.group}
    if args.limit:
        selected_rows = selected_rows[: args.limit]
    species_model = SpeciesModel(args.species_tree)
    iq_tree = Phylo.read(str(args.iqtree_species_tree), "newick")
    iq_species = {tip.name for tip in iq_tree.get_terminals()}
    rooted_species = set(species_model.by_species)
    if iq_species != rooted_species:
        raise ValueError(f"Rooted/IQ-TREE species mismatch: rooted_only={rooted_species-iq_species}, iq_only={iq_species-rooted_species}")

    family_rows = []
    duplication_rows = []
    loss_rows = []
    filtered_rows = []
    for selected in selected_rows:
        orthogroup = selected["orthogroup"]
        if orthogroup not in tree_rows:
            raise ValueError(f"No Step 10 tree summary row for {args.group} {orthogroup}")
        family, duplications, losses, filtered = reconcile_family(
            args.group, selected, tree_rows[orthogroup], species_model, args.output_dir,
            args.duplication_cost, args.loss_cost,
        )
        family_rows.append(family)
        duplication_rows.extend(duplications)
        loss_rows.extend(losses)
        filtered_rows.extend(filtered)

    family_fields = list(family_rows[0]) if family_rows else ["group", "orthogroup"]
    duplication_fields = list(duplication_rows[0]) if duplication_rows else [
        "group", "orthogroup", "gene_node_id", "species_node", "support_class", "optimal_root_robustness"
    ]
    loss_fields = list(loss_rows[0]) if loss_rows else ["group", "orthogroup", "loss_event_id", "lost_species_node"]
    write_tsv(args.output_dir / "family_reconciliation_summary.tsv", family_fields, family_rows)
    write_tsv(args.output_dir / "duplication_events.tsv", duplication_fields, duplication_rows)
    write_tsv(args.output_dir / "loss_events.tsv", loss_fields, loss_rows)
    write_tsv(
        args.output_dir / "alignment_filtered_sequences.tsv",
        ["group", "orthogroup", "tree_label", "species", "gene_id", "reason"],
        filtered_rows,
    )
    write_tsv(args.output_dir / "species_tree_branches.tsv", [
        "group", "species_node", "parent_species_node", "branch_type", "descendant_species_count", "descendant_species"
    ], species_model.branch_rows(args.group))

    dup_by_branch = defaultdict(list)
    loss_by_branch = defaultdict(list)
    for row in duplication_rows:
        dup_by_branch[row["species_node"]].append(row)
    for row in loss_rows:
        loss_by_branch[row["lost_species_node"]].append(row)
    branch_rows = []
    for branch in species_model.branch_rows(args.group):
        node = branch["species_node"]
        dups = dup_by_branch[node]
        losses = loss_by_branch[node]
        branch_rows.append({
            **branch,
            "families_with_duplication": len({row["orthogroup"] for row in dups}),
            "duplication_events": len(dups),
            "strong_support_duplications": sum(row["support_class"] == "strong" for row in dups),
            "root_robust_duplications": sum(number(row["optimal_root_robustness"]) == 1.0 for row in dups),
            "families_with_loss": len({row["orthogroup"] for row in losses}),
            "loss_events": len(losses),
        })
    write_tsv(args.output_dir / "branch_event_summary.tsv", list(branch_rows[0]), branch_rows)
    high_conf = [
        row for row in duplication_rows
        if row["support_class"] == "strong" and number(row["optimal_root_robustness"]) == 1.0
    ]
    write_tsv(args.output_dir / "high_confidence_duplication_events.tsv", duplication_fields, high_conf)
    summary = [{
        "group": args.group,
        "families": len(family_rows),
        "duplication_events": len(duplication_rows),
        "strong_support_duplications": sum(row["support_class"] == "strong" for row in duplication_rows),
        "root_robust_duplications": sum(number(row["optimal_root_robustness"]) == 1.0 for row in duplication_rows),
        "high_confidence_duplications": len(high_conf),
        "loss_events": len(loss_rows),
        "alignment_filtered_sequences": len(filtered_rows),
        "reinserted_iqtree_identical_sequences": sum(integer(row["reinserted_iqtree_identical_sequences"]) for row in family_rows),
        "families_with_multiple_optimal_roots": sum(integer(row["equally_optimal_roots"]) > 1 for row in family_rows),
        "median_equally_optimal_roots": median_int([integer(row["equally_optimal_roots"]) for row in family_rows]),
    }]
    write_tsv(args.output_dir / "step11_group_summary.tsv", list(summary[0]), summary)
    methods = f"""Step 11 gene-tree/species-tree reconciliation ({args.group})

Method: duplication-loss parsimony using classical LCA mapping. The rooted OrthoFinder/STRIDE species tree
provides the root; its taxon set was checked against the corresponding IQ-TREE species tree. For each unrooted
candidate gene tree, every edge was evaluated as a possible root. The selected root minimizes
{args.duplication_cost:g}*duplications + {args.loss_cost:g}*losses. All equally optimal roots are retained in the
root-candidate tables and used to calculate event rooting robustness.

Internal gene-tree support is read as SH-aLRT/UFBoot. A duplication is labelled strong only when SH-aLRT >=80
and UFBoot >=95. High-confidence events additionally require optimal-root robustness=1. This analysis models
duplication and loss only; it does not model incomplete lineage sorting, horizontal transfer, gene conversion,
assembly/annotation error, or uncertainty in the species tree. Event counts are hypotheses for downstream
interpretation, not direct observations or absolute dates.

Rooted species tree: {args.species_tree}
IQ-TREE species tree used for taxon verification: {args.iqtree_species_tree}
Step 10 input: {step10}
"""
    (args.output_dir / "METHODS_AND_LIMITATIONS.txt").write_text(methods)
    (args.output_dir / "COMPLETE").write_text(f"Step 11 {args.group} complete\n")


def median_int(values):
    if not values:
        return 0
    values = sorted(values)
    middle = len(values) // 2
    if len(values) % 2:
        return values[middle]
    return (values[middle - 1] + values[middle]) / 2


def concatenate_tables(paths, output):
    rows = []
    fields = None
    for path in paths:
        current = read_tsv(path)
        if current and fields is None:
            fields = list(current[0])
        rows.extend(current)
    write_tsv(output, fields or ["group", "orthogroup"], rows)
    return rows


def cmd_combine(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    group_dirs = {"bivalvia": args.bivalvia.resolve(), "gastropoda": args.gastropoda.resolve()}
    for group, directory in group_dirs.items():
        if not (directory / "COMPLETE").exists():
            raise FileNotFoundError(f"Incomplete Step 11 group result: {group} {directory}")
    family_rows = concatenate_tables(
        [directory / "family_reconciliation_summary.tsv" for directory in group_dirs.values()],
        args.output_dir / "family_reconciliation_summary.tsv",
    )
    duplication_rows = concatenate_tables(
        [directory / "duplication_events.tsv" for directory in group_dirs.values()],
        args.output_dir / "duplication_events.tsv",
    )
    loss_rows = concatenate_tables(
        [directory / "loss_events.tsv" for directory in group_dirs.values()],
        args.output_dir / "loss_events.tsv",
    )
    filtered_rows = concatenate_tables(
        [directory / "alignment_filtered_sequences.tsv" for directory in group_dirs.values()],
        args.output_dir / "alignment_filtered_sequences.tsv",
    )
    high_conf = concatenate_tables(
        [directory / "high_confidence_duplication_events.tsv" for directory in group_dirs.values()],
        args.output_dir / "high_confidence_duplication_events.tsv",
    )
    concatenate_tables(
        [directory / "branch_event_summary.tsv" for directory in group_dirs.values()],
        args.output_dir / "branch_event_summary.tsv",
    )
    concatenate_tables(
        [directory / "species_tree_branches.tsv" for directory in group_dirs.values()],
        args.output_dir / "species_tree_branches.tsv",
    )
    summary_rows = []
    for group in ("bivalvia", "gastropoda"):
        families = [row for row in family_rows if row["group"] == group]
        dups = [row for row in duplication_rows if row["group"] == group]
        losses = [row for row in loss_rows if row["group"] == group]
        filtered = [row for row in filtered_rows if row["group"] == group]
        highs = [row for row in high_conf if row["group"] == group]
        summary_rows.append({
            "group": group,
            "families": len(families),
            "duplication_events": len(dups),
            "strong_support_duplications": sum(row["support_class"] == "strong" for row in dups),
            "root_robust_duplications": sum(number(row["optimal_root_robustness"]) == 1.0 for row in dups),
            "high_confidence_duplications": len(highs),
            "loss_events": len(losses),
            "alignment_filtered_sequences": len(filtered),
            "reinserted_iqtree_identical_sequences": sum(integer(row["reinserted_iqtree_identical_sequences"]) for row in families),
            "families_with_multiple_optimal_roots": sum(integer(row["equally_optimal_roots"]) > 1 for row in families),
        })
    write_tsv(args.output_dir / "step11_summary.tsv", list(summary_rows[0]), summary_rows)
    priority_rows = []
    high_by_family = Counter((row["group"], row["orthogroup"]) for row in high_conf)
    strong_by_family = Counter((row["group"], row["orthogroup"]) for row in duplication_rows if row["support_class"] == "strong")
    for row in family_rows:
        key = (row["group"], row["orthogroup"])
        priority_rows.append({
            **row,
            "high_confidence_duplications": high_by_family[key],
            "strong_support_duplications_recount": strong_by_family[key],
            "reconciliation_priority_score": (
                high_by_family[key] * 100
                + strong_by_family[key] * 10
                + integer(row["duplications"])
                - integer(row["low_or_unknown_support_duplications"])
                - max(0, integer(row["equally_optimal_roots"]) - 1)
            ),
        })
    priority_rows.sort(key=lambda row: (-number(row["reconciliation_priority_score"]), row["group"], row["orthogroup"]))
    write_tsv(args.output_dir / "family_reconciliation_priorities.tsv", list(priority_rows[0]), priority_rows)
    provenance = [{"group": group, "result_directory": str(directory)} for group, directory in group_dirs.items()]
    write_tsv(args.output_dir / "input_result_directories.tsv", ["group", "result_directory"], provenance)
    methods = """Step 11 combined summary

This directory combines independent Bivalvia and Gastropoda duplication-loss LCA reconciliations. High-confidence
duplications require both strong IQ-TREE support (SH-aLRT >=80 and UFBoot >=95) and presence under every equally
optimal gene-tree root. Family priorities rank reconciliation evidence; they do not replace biological review.
See each group METHODS_AND_LIMITATIONS.txt for the full model and limitations.
"""
    (args.output_dir / "METHODS_AND_LIMITATIONS.txt").write_text(methods)
    manifest = []
    for path in sorted(args.output_dir.rglob("*")):
        if path.is_file() and path.name != "result_manifest.tsv":
            manifest.append({"relative_path": str(path.relative_to(args.output_dir)), "bytes": path.stat().st_size})
    write_tsv(args.output_dir / "result_manifest.tsv", ["relative_path", "bytes"], manifest)
    (args.output_dir / "COMPLETE").write_text("Step 11 combined reconciliation complete\n")


def build_parser():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    group = sub.add_parser("reconcile-group")
    group.add_argument("--group", choices=["bivalvia", "gastropoda"], required=True)
    group.add_argument("--step10", type=Path, required=True)
    group.add_argument("--species-tree", type=Path, required=True)
    group.add_argument("--iqtree-species-tree", type=Path, required=True)
    group.add_argument("--output-dir", type=Path, required=True)
    group.add_argument("--duplication-cost", type=float, default=1.0)
    group.add_argument("--loss-cost", type=float, default=1.0)
    group.add_argument("--limit", type=int, default=0)
    group.set_defaults(func=cmd_reconcile_group)
    combine = sub.add_parser("combine")
    combine.add_argument("--bivalvia", type=Path, required=True)
    combine.add_argument("--gastropoda", type=Path, required=True)
    combine.add_argument("--output-dir", type=Path, required=True)
    combine.set_defaults(func=cmd_combine)
    return parser


if __name__ == "__main__":
    arguments = build_parser().parse_args()
    arguments.func(arguments)
