"""Graph schema inferred from biological network JSON (single source of truth)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

# BiologicalGraph imported lazily in functions to avoid circular imports.

# Types that are never directly measured in omics tables.
DEFAULT_NON_OBSERVABLE = frozenset({"R", "P"})

# Graph node type -> expected wide CSV filename (stable on disk).
CSV_NAME_OVERRIDES: Dict[str, str] = {
    "miRNA": "mirna.csv",
    "Enhancer": "enhancer.csv",
    "Promoter": "promoter.csv",
    "Promoter/Enhancer": "promoter_enhancer.csv",
    "Promoter_Enhancer": "promoter_enhancer.csv",
    "protein_complex": "protein_complex.csv",
    "snp": "snp.csv",
    "gene": "gene.csv",
    "transcript": "transcript.csv",
    "protein": "protein.csv",
    "mirna": "mirna.csv",
    "enhancer": "enhancer.csv",
    "promoter": "promoter.csv",
    "promoter_enhancer": "promoter_enhancer.csv",
}


def csv_for_omic_type(omic_type: str) -> str:
    """Map store / graph omic label to on-disk wide CSV filename."""
    return CSV_NAME_OVERRIDES.get(omic_type, f"{omic_type.lower().replace('/', '_')}.csv")


# Node types that receive entity suffixes when mapping tabular features to graph IDs.
SUFFIX_NODE_TYPES: Dict[str, str] = {
    "transcript": "_transcript",
    "protein": "_protein",
}


@dataclass
class GraphSchema:
    """Node/edge vocabulary extracted from a biological graph."""

    node_types: List[str]
    edge_types: List[Tuple[str, str]]
    relation_types: List[str]
    non_observable_types: List[str] = field(default_factory=lambda: sorted(DEFAULT_NON_OBSERVABLE))
    omics_file_map: Dict[str, str] = field(default_factory=dict)

    @property
    def observable_types(self) -> List[str]:
        non_obs = set(self.non_observable_types)
        return [t for t in self.node_types if t not in non_obs]

    def csv_for_type(self, node_type: str) -> str:
        return CSV_NAME_OVERRIDES.get(node_type, f"{node_type.lower().replace('/', '_')}.csv")

    def build_omics_file_map(self, included_types: Iterable[str] | None = None) -> Dict[str, str]:
        types = list(included_types) if included_types is not None else self.observable_types
        return {self.csv_for_type(t): t for t in types if t in self.node_types}


def infer_node_types(graph) -> List[str]:
    types: Set[str] = set()
    for meta in graph.vertices.values():
        otype = meta.get("omic_type", "unknown")
        if otype:
            types.add(str(otype))
    return sorted(types)


def infer_edge_types(graph) -> List[Tuple[str, str]]:
    pairs: Set[Tuple[str, str]] = set()
    for edge in graph.edges:
        s, t = edge["start_vertex"], edge["end_vertex"]
        if s not in graph.vertices or t not in graph.vertices:
            continue
        st = str(graph.vertices[s].get("omic_type", "unknown"))
        tt = str(graph.vertices[t].get("omic_type", "unknown"))
        pairs.add((st, tt))
    return sorted(pairs)


def infer_relation_types(graph) -> List[str]:
    rels: Set[str] = set()
    for edge in graph.edges:
        int_info = edge.get("int_info") or {}
        rel = int_info.get("type")
        if rel:
            rels.add(str(rel))
    return sorted(rels)


def schema_from_graph(graph) -> GraphSchema:
    node_types = infer_node_types(graph)
    edge_types = infer_edge_types(graph)
    relation_types = infer_relation_types(graph)
    non_obs = sorted(t for t in DEFAULT_NON_OBSERVABLE if t in node_types)
    omics_map = {CSV_NAME_OVERRIDES.get(t, f"{t.lower().replace('/', '_')}.csv"): t
                 for t in node_types if t not in non_obs}
    return GraphSchema(
        node_types=node_types,
        edge_types=edge_types,
        relation_types=relation_types,
        non_observable_types=non_obs,
        omics_file_map=omics_map,
    )


def load_graph_schema(path: str | Path) -> GraphSchema:
    from universal_graph_mapper.biological_graph import BiologicalGraph
    graph = BiologicalGraph.from_json_path(path)
    return schema_from_graph(graph)


def store_config_slug(
    omics_inclusion: Iterable[str],
    *,
    prune_k_hops: int = 0,
    prune_direction: str = "out",
    exclude_node_types: Iterable[str] | None = None,
) -> str:
    """Deterministic store name from upstream preprocessing config."""
    parts = sorted(str(t).replace("/", "_").replace(" ", "_").lower() for t in omics_inclusion)
    slug = "_".join(parts) if parts else "all_omics"
    if prune_k_hops > 0:
        slug += f"_prune{prune_k_hops}{prune_direction[0]}"
    excluded = sorted(exclude_node_types or [])
    if excluded:
        slug += "_excl_" + "_".join(t.replace("/", "_").lower() for t in excluded)
    return slug
