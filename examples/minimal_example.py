#!/usr/bin/env python3
"""Minimal API example: map a list of feature IDs to universal graph nodes.

Uses only packaged data under data/ (see config/paths.yaml).
"""

from pathlib import Path

from universal_graph_mapper.config import load_paths
from universal_graph_mapper.biological_graph import BiologicalGraph
from universal_graph_mapper.json_loader import vertex_ids_by_type
from universal_graph_mapper.mapping.pipeline import MappingStrategy, feature_rows_to_frame, run_legacy_then_optimize
from universal_graph_mapper.mapping.resolver import OptimizedResolver
from universal_graph_mapper.mapping.sources import make_sources_for_profile

PKG_ROOT = Path(__file__).resolve().parent.parent
paths = load_paths(PKG_ROOT / "config" / "paths.yaml")
paths.validate()

print("Loading graph:", paths.universal_graph)
graph = BiologicalGraph.from_json_path(paths.universal_graph)
gene_entities = set(vertex_ids_by_type(graph).get("gene", []))

features = ["ENSG00000141510", "TP53", "P04637", "Q07973"]
resolver = OptimizedResolver(sources=make_sources_for_profile("overall"), strategy_label="overall")

rows = run_legacy_then_optimize(
    features,
    node_type="gene",
    graph_entities=gene_entities,
    resolver=resolver,
    strategy=MappingStrategy.OPTIMIZED,
)
print(feature_rows_to_frame(rows).to_string(index=False))
