"""Prepare feature_map with optional optimized recovery (legacy preserved)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, Set

import pandas as pd

from universal_graph_mapper.graph_schema import csv_for_omic_type
from universal_graph_mapper.constants import ID_COLUMNS
from universal_graph_mapper.biological_graph import BiologicalGraph
from universal_graph_mapper.json_loader import vertex_ids_by_type
from universal_graph_mapper.feature_mapping import apply_entity_suffix
from universal_graph_mapper.mapping.pipeline import (
    MappingStrategy,
    feature_rows_to_frame,
    run_legacy_then_optimize,
)
from universal_graph_mapper.mapping.resolver import OptimizedResolver
from universal_graph_mapper.mapping.sources import make_sources_for_profile


def _is_id_column(col: str) -> bool:
    return col.strip().lower() in ID_COLUMNS


def _entities_for_type(graph: BiologicalGraph, node_type: str) -> Set[str]:
    by_type = vertex_ids_by_type(graph)
    return set(by_type.get(node_type, []))


def prepare_multi_cohort_raw_with_strategy(
    raw_dir: Path,
    omics_files: Dict[str, str],
    *,
    graph_path: Path | None = None,
    omics_read_dir: Path | None = None,
    omics_type_overrides: Dict[str, str] | None = None,
    mapping_strategy: str | MappingStrategy = MappingStrategy.OPTIMIZED,
    use_optimized_mapping: bool | None = None,
    mapping_profile: str = "overall",
    write_optimized_sidecar: bool = True,
) -> dict:
    """Build feature_map.csv using legacy and/or optimized mapping.

    Parameters
    ----------
    mapping_strategy:
        ``legacy`` | ``optimized`` | ``compare``.
        ``optimized`` (default) = legacy first, then recover unmapped via ``mapping_profile``.
        ``legacy`` = original exact-match behaviour only.
        ``compare`` = same as optimized for the written feature_map, but also
        writes separate legacy / optimized frames under ``results/mapping/``.
    mapping_profile:
        ``biomart`` | ``hgnc`` | ``uniprot`` | ``overall`` (default overall =
        BioMart + HGNC + UniProt).
    use_optimized_mapping:
        Convenience flag; if False forces legacy, if True forces optimized.
        Overrides ``mapping_strategy`` when not None.
    """
    if use_optimized_mapping is False:
        strategy = MappingStrategy.LEGACY
    elif use_optimized_mapping is True:
        strategy = MappingStrategy.OPTIMIZED
    else:
        strategy = MappingStrategy(mapping_strategy)

    raw_dir = Path(raw_dir)
    read_dir = Path(omics_read_dir) if omics_read_dir else raw_dir
    if graph_path is None:
        from universal_graph_mapper.config import load_paths
        graph_path = load_paths().universal_graph
    graph_path = Path(graph_path)
    graph = BiologicalGraph.from_json_path(graph_path)
    overrides = omics_type_overrides or {}
    resolver = OptimizedResolver(
        sources=make_sources_for_profile(mapping_profile),
        strategy_label=mapping_profile,
    )

    all_rows = []
    measured: Set[str] = set()
    graph_entities_by_type: dict[str, set[str]] = {}

    for fname, graph_type in omics_files.items():
        map_type = overrides.get(fname, graph_type)
        entities = _entities_for_type(graph, map_type)
        graph_entities_by_type[map_type] = entities
        path = read_dir / fname
        if not path.exists():
            raise FileNotFoundError(f"Missing {path}")
        header = [c for c in pd.read_csv(path, nrows=0).columns if not _is_id_column(c)]
        rows = run_legacy_then_optimize(
            header,
            node_type=map_type,
            graph_entities=entities,
            omics_file=fname,
            resolver=resolver,
            strategy=strategy if strategy is not MappingStrategy.COMPARE else MappingStrategy.OPTIMIZED,
        )
        all_rows.extend(rows)
        for r in rows:
            if r.entity_id:
                measured.add(r.entity_id)
            for eid in r.entity_ids:
                measured.add(eid)

    fmap = feature_rows_to_frame(all_rows)
    # Canonical columns expected by the rest of the pipeline.
    out = fmap[
        ["omics_file", "feature_col", "entity_id", "node_type", "in_graph"]
    ].copy()
    out.to_csv(raw_dir / "feature_map.csv", index=False)

    if write_optimized_sidecar:
        fmap.to_csv(raw_dir / "feature_map_detailed.csv", index=False)

    with open(raw_dir / "measured_entities.json", "w") as f:
        json.dump(sorted(measured), f, indent=2)

    n_in = int(out["in_graph"].sum())
    print(f"  feature_map.csv: {len(out)} rows ({n_in} in-graph) [strategy={strategy.value}]")
    print(f"  measured_entities.json: {len(measured)} entities")
    return {
        "n_features": len(out),
        "n_in_graph": n_in,
        "n_measured": len(measured),
        "strategy": strategy.value,
        "detailed": fmap,
        "graph_entities_by_type": {k: len(v) for k, v in graph_entities_by_type.items()},
    }


def default_omics_files(included_types: Iterable[str]) -> Dict[str, str]:
    return {csv_for_omic_type(t): t for t in included_types}


# Re-export suffix helper so callers needing it stay on one import path.
__all__ = [
    "prepare_multi_cohort_raw_with_strategy",
    "default_omics_files",
    "apply_entity_suffix",
]
