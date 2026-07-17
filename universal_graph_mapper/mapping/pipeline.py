"""Legacy-first then optimized recovery for gene/transcript/protein features."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from universal_graph_mapper.feature_mapping import apply_entity_suffix
from universal_graph_mapper.mapping.id_detect import detect_identifier_type
from universal_graph_mapper.mapping.resolver import (
    OPTIMIZABLE_NODE_TYPES,
    MappingResult,
    OptimizedResolver,
)
from universal_graph_mapper.mapping.sources import MappingSource


class MappingStrategy(str, Enum):
    LEGACY = "legacy"
    OPTIMIZED = "optimized"
    COMPARE = "compare"  # run both; prefer legacy when present, recover with optimized


@dataclass
class FeatureMapRow:
    omics_file: str
    feature_col: str
    node_type: str
    entity_id: str
    entity_ids: list[str] = field(default_factory=list)
    in_graph: bool = False
    strategy: str = "legacy"
    detected_type: str = ""
    mapping_path: str = ""
    uniprot_ids: list[str] = field(default_factory=list)
    status: str = ""
    failure_reason: str = ""


def legacy_map_column(col: str, node_type: str, entities: set[str]) -> tuple[str | None, bool]:
    """Pure-function clone of prepare_multi_cohort_raw legacy lookup (untouched source)."""
    candidates = [col, col.strip(), apply_entity_suffix(col, node_type)]
    if node_type == "gene":
        for c in list(candidates):
            if c.endswith("_transcript"):
                candidates.append(c[: -len("_transcript")])
    seen: set[str] = set()
    for cand in candidates:
        if cand in seen:
            continue
        seen.add(cand)
        if cand in entities:
            return cand, True
    return None, False


def run_legacy_then_optimize(
    features: Sequence[str],
    *,
    node_type: str,
    graph_entities: set[str],
    omics_file: str = "",
    resolver: OptimizedResolver | None = None,
    strategy: MappingStrategy | str = MappingStrategy.OPTIMIZED,
) -> list[FeatureMapRow]:
    """Map features with legacy first; optimize only unmapped (never overwrite)."""
    strat = MappingStrategy(strategy)
    resolver = resolver or OptimizedResolver()
    rows: list[FeatureMapRow] = []

    for col in features:
        detected = detect_identifier_type(col)
        legacy_eid, legacy_ok = legacy_map_column(col, node_type, graph_entities)

        if node_type not in OPTIMIZABLE_NODE_TYPES:
            rows.append(
                FeatureMapRow(
                    omics_file=omics_file,
                    feature_col=col,
                    node_type=node_type,
                    entity_id=legacy_eid or "",
                    entity_ids=[legacy_eid] if legacy_eid else [],
                    in_graph=legacy_ok,
                    strategy="legacy",
                    detected_type=detected.value,
                    mapping_path="Legacy mapper" if legacy_ok else "Not mapped",
                    status="success" if legacy_ok else "unmapped",
                    failure_reason="" if legacy_ok else "out_of_scope_or_legacy_miss",
                )
            )
            continue

        if legacy_ok and legacy_eid:
            # Never overwrite a valid legacy UniProt/graph mapping.
            rows.append(
                FeatureMapRow(
                    omics_file=omics_file,
                    feature_col=col,
                    node_type=node_type,
                    entity_id=legacy_eid,
                    entity_ids=[legacy_eid],
                    in_graph=True,
                    strategy="legacy",
                    detected_type=detected.value,
                    mapping_path="Legacy mapper",
                    uniprot_ids=[_uniprot_from_entity(legacy_eid)],
                    status="success",
                )
            )
            if strat is MappingStrategy.LEGACY:
                continue
            # OPTIMIZED / COMPARE also keep legacy hit.
            continue

        if strat is MappingStrategy.LEGACY:
            rows.append(
                FeatureMapRow(
                    omics_file=omics_file,
                    feature_col=col,
                    node_type=node_type,
                    entity_id="",
                    in_graph=False,
                    strategy="legacy",
                    detected_type=detected.value,
                    mapping_path="Not mapped",
                    status="unmapped",
                    failure_reason="legacy_miss",
                )
            )
            continue

        # Step 3: optimized recovery only on unmapped features.
        opt = resolver.resolve_to_graph_entities(col, node_type, graph_entities)
        primary = opt.entity_ids[0] if opt.entity_ids else ""
        # Prefer resolver source profile name when set (biomart/hgnc/uniprot/overall).
        strat_label = getattr(resolver, "strategy_label", None) or "optimized"
        rows.append(
            FeatureMapRow(
                omics_file=omics_file,
                feature_col=col,
                node_type=node_type,
                entity_id=primary,
                entity_ids=list(opt.entity_ids),
                in_graph=bool(opt.entity_ids),
                strategy=strat_label,
                detected_type=opt.detected_type.value,
                mapping_path=opt.mapping_path,
                uniprot_ids=list(opt.uniprot_ids),
                status=opt.status if opt.entity_ids else "unmapped",
                failure_reason=opt.failure_reason,
            )
        )
    return rows


def map_features_for_dataset(
    omics_columns: dict[str, Sequence[str]],
    *,
    omics_node_types: dict[str, str],
    graph_entities_by_type: dict[str, set[str]],
    strategy: MappingStrategy | str = MappingStrategy.OPTIMIZED,
    sources: Sequence[MappingSource] | None = None,
) -> pd.DataFrame:
    """Map every omics file/column set; returns a feature_map-like DataFrame."""
    resolver = OptimizedResolver(sources=sources) if sources is not None else OptimizedResolver()
    all_rows: list[FeatureMapRow] = []
    for fname, cols in omics_columns.items():
        node_type = omics_node_types[fname]
        entities = graph_entities_by_type.get(node_type, set())
        all_rows.extend(
            run_legacy_then_optimize(
                cols,
                node_type=node_type,
                graph_entities=entities,
                omics_file=fname,
                resolver=resolver,
                strategy=strategy,
            )
        )
    return feature_rows_to_frame(all_rows)


def feature_rows_to_frame(rows: Iterable[FeatureMapRow]) -> pd.DataFrame:
    records = []
    for r in rows:
        records.append(
            {
                "omics_file": r.omics_file,
                "feature_col": r.feature_col,
                "entity_id": r.entity_id,
                "entity_ids": "|".join(r.entity_ids),
                "node_type": r.node_type,
                "in_graph": r.in_graph,
                "strategy": r.strategy,
                "detected_type": r.detected_type,
                "mapping_path": r.mapping_path,
                "uniprot_ids": "|".join(r.uniprot_ids),
                "status": r.status,
                "failure_reason": r.failure_reason,
            }
        )
    return pd.DataFrame.from_records(records)


def _uniprot_from_entity(entity_id: str) -> str:
    eid = str(entity_id)
    for suffix in ("_transcript", "_protein"):
        if eid.endswith(suffix):
            return eid[: -len(suffix)]
    return eid
