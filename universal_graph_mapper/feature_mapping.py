"""Feature column → graph entity mapping with transcript/protein suffix rules."""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple

from universal_graph_mapper.graph_schema import SUFFIX_NODE_TYPES

# Graph entity IDs for transcript/protein carry these suffixes; gene IDs do not.
_CROSS_OMICS_SUFFIXES = tuple(SUFFIX_NODE_TYPES.values())  # ("_transcript", "_protein")


def apply_entity_suffix(col: str, node_type: str) -> str:
    """Only transcript and protein omics receive graph ID suffixes."""
    suffix = SUFFIX_NODE_TYPES.get(node_type)
    if suffix and not col.endswith(suffix):
        return f"{col}{suffix}"
    return col


def _unsuffix_candidates(entity_id: str) -> List[str]:
    """Yield gene-style base IDs from suffixed transcript/protein entity IDs."""
    out = [entity_id]
    for suf in _CROSS_OMICS_SUFFIXES:
        if entity_id.endswith(suf):
            out.append(entity_id[: -len(suf)])
    return out


def resolve_entity_id(
    col: str,
    node_type: str,
    entity_map: Dict[str, int],
    type_entity_map: Dict[str, int],
    global_feature_map: Dict[str, str],
    file_feature_map: Dict[str, str],
) -> Tuple[Optional[str], Optional[int]]:
    """Map a feature column name to a registry entity_id and node index."""
    candidates: List[str] = []
    if col in file_feature_map:
        candidates.append(file_feature_map[col])
    if col in global_feature_map:
        candidates.append(global_feature_map[col])
    candidates.extend([col, col.strip(), apply_entity_suffix(col, node_type)])

    # transcript.csv → gene override: feature_map entities are ``*_transcript`` but
    # the pruned gene registry uses unsuffixed gene IDs (same base).
    expanded: List[str] = []
    for cand in candidates:
        expanded.extend(_unsuffix_candidates(cand))
        if node_type in SUFFIX_NODE_TYPES:
            expanded.append(apply_entity_suffix(cand, node_type))
    candidates = expanded

    seen = set()
    ordered: List[str] = []
    for cand in candidates:
        if cand in seen:
            continue
        seen.add(cand)
        ordered.append(cand)

    # Prefer type-filtered registry matches (critical for transcript.csv → gene).
    for cand in ordered:
        if cand in type_entity_map:
            return cand, type_entity_map[cand]
    for cand in ordered:
        if cand in entity_map:
            return cand, entity_map[cand]
    return None, None


def adapt_measured_entities_for_overrides(
    entities: Iterable[str],
    *,
    selected_omics: Iterable[str] | None = None,
    omics_type_overrides: Dict[str, str] | None = None,
) -> List[str]:
    """Expand measured entity seeds so gene selections prune from gene IDs.

    When ``transcript.csv`` is remapped to ``gene``, feature_map seeds are
    ``*_transcript``; k-hop pruning with ``omics_inclusion`` containing ``gene``
    needs the unsuffixed gene IDs as seeds.
    """
    selected = set(selected_omics or [])
    overrides = omics_type_overrides or {}
    needs_gene_bases = "gene" in selected or "gene" in overrides.values()
    out: set[str] = set()
    for ent in entities:
        e = str(ent).strip()
        if not e or e.lower() == "nan":
            continue
        out.add(e)
        if needs_gene_bases:
            for base in _unsuffix_candidates(e):
                out.add(base)
    return sorted(out)
