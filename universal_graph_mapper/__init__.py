"""Map omics feature IDs onto a universal biological graph (MOGI-compatible)."""

from universal_graph_mapper.config import MapperPaths, load_paths
from universal_graph_mapper.mapping import (
    MappingStrategy,
    map_features_for_dataset,
    prepare_multi_cohort_raw_with_strategy,
    run_legacy_then_optimize,
)

__version__ = "1.0.0"
__all__ = [
    "MapperPaths",
    "load_paths",
    "MappingStrategy",
    "map_features_for_dataset",
    "prepare_multi_cohort_raw_with_strategy",
    "run_legacy_then_optimize",
]
