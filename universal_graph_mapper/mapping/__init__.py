"""Feature → universal graph mapping (legacy + optimized)."""

from universal_graph_mapper.mapping.id_detect import (
    IdentifierType,
    detect_identifier_type,
    normalize_identifier,
)
from universal_graph_mapper.mapping.pipeline import (
    MappingStrategy,
    feature_rows_to_frame,
    map_features_for_dataset,
    run_legacy_then_optimize,
)
from universal_graph_mapper.mapping.prepare import prepare_multi_cohort_raw_with_strategy
from universal_graph_mapper.mapping.resolver import MappingResult, OptimizedResolver
from universal_graph_mapper.mapping.sources import (
    BioMartSource,
    HgncSource,
    UniProtIdMappingSource,
    make_source,
    make_sources_for_profile,
)

__all__ = [
    "IdentifierType",
    "detect_identifier_type",
    "normalize_identifier",
    "MappingStrategy",
    "MappingResult",
    "OptimizedResolver",
    "feature_rows_to_frame",
    "map_features_for_dataset",
    "run_legacy_then_optimize",
    "prepare_multi_cohort_raw_with_strategy",
    "BioMartSource",
    "HgncSource",
    "UniProtIdMappingSource",
    "make_source",
    "make_sources_for_profile",
]
