"""Biological graph JSON schema (MOGI-compatible backbone format)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Set

GRAPH_JSON_SCHEMA_VERSION = "1.0"


@dataclass
class GraphVertex:
    """One node in a universal biological graph."""

    node_id: str
    omic_type: str
    label: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class GraphEdge:
    """Directed edge between two vertices."""

    start_vertex: str
    end_vertex: str
    relation_type: str = "interacts"
    metadata: Dict[str, Any] = field(default_factory=dict)


def validate_graph_json(data: Dict[str, Any]) -> None:
    """Raise ValueError if required graph JSON keys are missing."""
    if "vertices" not in data or "edges" not in data:
        raise ValueError("Graph JSON must contain 'vertices' and 'edges'")
    if not isinstance(data["vertices"], dict):
        raise ValueError("'vertices' must be a dict")
    if not isinstance(data["edges"], list):
        raise ValueError("'edges' must be a list")

    vertex_ids: Set[str] = set(data["vertices"].keys())
    for i, edge in enumerate(data["edges"]):
        for key in ("start_vertex", "end_vertex"):
            if key not in edge:
                raise ValueError(f"Edge {i} missing '{key}'")
            if edge[key] not in vertex_ids:
                raise ValueError(f"Edge {i} references unknown vertex '{edge[key]}'")


def normalize_omic_type(raw: str) -> str:
    """Map external / MOGI omic labels to canonical lowercase types (legacy helper)."""
    aliases = {
        "mirna": "mirna",
        "microRNA": "mirna",
        "miRNA": "mirna",
        "Enhancer": "enhancer",
        "Promoter": "promoter",
        "Promoter/Enhancer": "promoter_enhancer",
        "Promoter_Enhancer": "promoter_enhancer",
        "protein_complex": "protein_complex",
        "gene": "gene",
        "transcript": "transcript",
        "protein": "protein",
        "snp": "snp",
        "R": "R",
        "P": "P",
    }
    return aliases.get(raw, raw)


def infer_layer_order(node_types: list[str]) -> dict[str, int]:
    """Assign layout layers from graph-native node type names."""
    defaults = {
        "snp": 0, "enhancer": 0, "Enhancer": 0, "promoter": 0, "Promoter": 0,
        "promoter_enhancer": 0, "Promoter/Enhancer": 0, "mirna": 0, "miRNA": 0,
        "gene": 1, "transcript": 2, "protein": 3, "protein_complex": 3,
        "R": 4, "P": 5,
    }
    order = {}
    for t in node_types:
        order[t] = defaults.get(t, 2)
    return order
