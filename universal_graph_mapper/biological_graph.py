"""In-memory biological graph (MOGI-compatible JSON structure)."""

from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set

from universal_graph_mapper.graph_validation import validate_graph_json


@dataclass
class BiologicalGraph:
    """Universal biological graph as vertices + edges dicts."""

    vertices: Dict[str, Dict[str, Any]]
    edges: List[Dict[str, Any]]
    directed: bool = True
    schema_version: str = "1.0"

    def copy(self) -> "BiologicalGraph":
        return BiologicalGraph(
            vertices=deepcopy(self.vertices),
            edges=deepcopy(self.edges),
            directed=self.directed,
            schema_version=self.schema_version,
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "directed": self.directed,
            "vertices": self.vertices,
            "edges": self.edges,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "BiologicalGraph":
        validate_graph_json(data)
        return cls(
            vertices=dict(data["vertices"]),
            edges=list(data["edges"]),
            directed=bool(data.get("directed", True)),
            schema_version=data.get("schema_version", "1.0"),
        )

    @classmethod
    def from_json_path(cls, path: str | Path) -> "BiologicalGraph":
        with open(path) as f:
            return cls.from_dict(json.load(f))

    def save_json(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    def node_types(self) -> Set[str]:
        types: Set[str] = set()
        for v in self.vertices.values():
            otype = v.get("omic_type", "unknown")
            if otype:
                types.add(str(otype))
        return types

    def nodes_of_type(self, omic_type: str) -> List[str]:
        want = str(omic_type)
        return [
            nid for nid, v in self.vertices.items()
            if str(v.get("omic_type", "")) == want
        ]

    def measured_nodes_by_type(
        self,
        measured: Iterable[str],
    ) -> Dict[str, List[str]]:
        measured_set = set(measured)
        out: Dict[str, List[str]] = {}
        for nid in measured_set:
            if nid not in self.vertices:
                continue
            otype = str(self.vertices[nid].get("omic_type", "unknown"))
            out.setdefault(otype, []).append(nid)
        return out
