"""Graph JSON helpers (universal graph loader utilities)."""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Set

from universal_graph_mapper.biological_graph import BiologicalGraph


def load_biological_graph(path: str) -> BiologicalGraph:
    return BiologicalGraph.from_json_path(path)


def vertex_ids_by_type(
    graph: BiologicalGraph,
    *,
    exclude_types: Set[str] | None = None,
) -> Dict[str, List[str]]:
    """Group graph vertex IDs by omic_type."""
    if exclude_types is None:
        exclude_types = set()
    node_ids_by_type: Dict[str, List[str]] = defaultdict(list)
    for vid, meta in graph.vertices.items():
        ntype = str(meta.get("omic_type", "unknown"))
        if ntype in exclude_types:
            continue
        node_ids_by_type[ntype].append(vid)
    return dict(node_ids_by_type)
