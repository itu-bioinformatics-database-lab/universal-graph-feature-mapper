"""Multi-step resolver: feature ID → UniProt accessions with provenance."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Sequence

from universal_graph_mapper.graph_schema import SUFFIX_NODE_TYPES
from universal_graph_mapper.mapping.id_detect import (
    IdentifierType,
    detect_identifier_type,
    normalize_identifier,
)
from universal_graph_mapper.mapping.sources import BioMartSource, MappingSource

# Node types this mapper is allowed to optimize.
OPTIMIZABLE_NODE_TYPES = frozenset({"gene", "transcript", "protein"})

# Hops used to recover in-network UniProt siblings for out-of-network IDs.
_BRIDGE_HOP_TYPES = frozenset(
    {
        IdentifierType.ENSEMBL_GENE,
        IdentifierType.GENE_SYMBOL,
    }
)


@dataclass
class MappingResult:
    """Provenance-preserving mapping of one feature column."""

    feature_id: str
    detected_type: IdentifierType
    uniprot_ids: list[str] = field(default_factory=list)
    entity_ids: list[str] = field(default_factory=list)
    mapping_path: str = "Not mapped"
    strategy: str = "optimized"
    status: str = "unmapped"
    failure_reason: str = ""
    attempts: list[str] = field(default_factory=list)

    @property
    def mapped(self) -> bool:
        return self.status == "success" and bool(self.uniprot_ids)


@dataclass
class _Node:
    identifier: str
    id_type: IdentifierType
    path: tuple[str, ...]


class OptimizedResolver:
    """BFS over plugged-in mapping sources until UniProt accessions are reached."""

    def __init__(
        self,
        sources: Sequence[MappingSource] | None = None,
        *,
        max_depth: int = 4,
        strategy_label: str = "optimized",
    ):
        self.sources: list[MappingSource] = list(sources) if sources is not None else [BioMartSource()]
        self.max_depth = max_depth
        self.strategy_label = strategy_label

    def resolve_to_uniprot(self, feature_id: str) -> MappingResult:
        detected = detect_identifier_type(feature_id)
        norm = normalize_identifier(feature_id)
        result = MappingResult(
            feature_id=str(feature_id),
            detected_type=detected,
            strategy=self.strategy_label,
        )
        if detected in (IdentifierType.EMPTY, IdentifierType.OTHER):
            result.failure_reason = (
                "empty_identifier" if detected is IdentifierType.EMPTY else "unrecognized_identifier_type"
            )
            result.mapping_path = "Not mapped"
            return result

        if detected is IdentifierType.UNIPROT:
            uid = _strip_isoform(norm)
            result.uniprot_ids = [uid]
            result.mapping_path = "Already UniProt"
            result.status = "success"
            return result

        found: dict[str, tuple[str, ...]] = {}
        visited: set[tuple[str, IdentifierType]] = set()
        queue: deque[_Node] = deque([_Node(norm, detected, (detected.value,))])
        attempts: list[str] = []

        while queue:
            node = queue.popleft()
            key = (node.identifier, node.id_type)
            if key in visited:
                continue
            visited.add(key)
            if len(node.path) - 1 > self.max_depth:
                continue

            if node.id_type is IdentifierType.UNIPROT:
                uid = _strip_isoform(node.identifier)
                if uid not in found:
                    found[uid] = node.path
                continue

            expanded = False
            for source in self.sources:
                hits = source.lookup(node.identifier, node.id_type)
                if not hits:
                    continue
                expanded = True
                for tid, t_type in hits:
                    attempts.append(f"{source.name}:{node.id_type.value}->{t_type.value}:{tid}")
                    queue.append(_Node(tid, t_type, node.path + (t_type.value,)))
            if not expanded:
                attempts.append(f"dead_end:{node.id_type.value}:{node.identifier}")

        result.attempts = attempts[:50]
        if not found:
            result.failure_reason = "missing_mapping"
            result.mapping_path = "Not mapped"
            return result

        # Stable sorted UniProt list; path uses shortest discovery path summary.
        result.uniprot_ids = sorted(found.keys())
        best_path = min(found.values(), key=len)
        result.mapping_path = _format_path(best_path)
        result.status = "success"
        return result

    def resolve_to_graph_entities(
        self,
        feature_id: str,
        node_type: str,
        graph_entities: set[str],
    ) -> MappingResult:
        """Map feature → UniProt → graph entity IDs present in the skeleton.

        If the direct UniProt accession is not in the network, bridge through
        Ensembl gene / gene symbol to recover **sibling** UniProt IDs that are
        present (one-to-many safe)::

            UniProtB (data, not in graph)
              → ENSG1 / SYMBOL
              → {UniProtA, UniProtB, ...}
              → keep UniProtA if in graph
        """
        if node_type not in OPTIMIZABLE_NODE_TYPES:
            return MappingResult(
                feature_id=str(feature_id),
                detected_type=detect_identifier_type(feature_id),
                strategy=self.strategy_label,
                status="skipped",
                failure_reason="out_of_scope_node_type",
                mapping_path="Out of scope",
            )

        result = self.resolve_to_uniprot(feature_id)
        if not result.mapped:
            return result

        entity_ids = _match_uniprots_to_entities(result.uniprot_ids, node_type, graph_entities)
        result.entity_ids = entity_ids
        if entity_ids:
            result.status = "success"
            return result

        # Direct UniProt(s) missed the graph — try sibling recovery via ENSG/symbol.
        bridged_uids, bridge_path, bridge_attempts = self._bridge_to_network_uniprots(
            result.uniprot_ids,
            node_type=node_type,
            graph_entities=graph_entities,
        )
        result.attempts = (result.attempts + bridge_attempts)[:80]
        if bridged_uids:
            # Prefer in-network siblings; keep original out-of-network IDs in the list
            # for provenance (sorted unique).
            all_uids = sorted(set(result.uniprot_ids) | set(bridged_uids))
            result.uniprot_ids = all_uids
            result.entity_ids = _match_uniprots_to_entities(
                bridged_uids, node_type, graph_entities
            )
            result.mapping_path = bridge_path
            result.status = "success"
            result.failure_reason = ""
            return result

        result.status = "unmapped"
        result.failure_reason = "uniprot_not_in_network"
        result.mapping_path = f"{result.mapping_path} (not in network)"
        return result

    def _bridge_to_network_uniprots(
        self,
        seed_uniprots: Sequence[str],
        *,
        node_type: str,
        graph_entities: set[str],
    ) -> tuple[list[str], str, list[str]]:
        """UniProt → Ensembl/symbol → UniProt siblings that hit the network."""
        seeds = {_strip_isoform(u) for u in seed_uniprots}
        if not seeds or not self.sources:
            return [], "", []

        found: dict[str, tuple[str, ...]] = {}  # in-network UniProt → path
        visited: set[tuple[str, IdentifierType]] = set()
        attempts: list[str] = []
        queue: deque[_Node] = deque(
            _Node(uid, IdentifierType.UNIPROT, (IdentifierType.UNIPROT.value,)) for uid in sorted(seeds)
        )

        while queue:
            node = queue.popleft()
            key = (node.identifier, node.id_type)
            if key in visited:
                continue
            visited.add(key)
            depth = len(node.path) - 1
            if depth > self.max_depth:
                continue

            if node.id_type is IdentifierType.UNIPROT:
                uid = _strip_isoform(node.identifier)
                # Record only non-seed UniProts that sit on the graph (sibling recovery).
                if uid not in seeds and _match_uniprots_to_entities([uid], node_type, graph_entities):
                    if uid not in found or len(node.path) < len(found[uid]):
                        found[uid] = node.path
                # Expand seed (and intermediate) UniProts only into ENSG / symbol bridges.
                if depth >= self.max_depth:
                    continue
                for source in self.sources:
                    for tid, t_type in source.lookup(uid, IdentifierType.UNIPROT):
                        if t_type not in _BRIDGE_HOP_TYPES:
                            continue
                        attempts.append(
                            f"bridge:{source.name}:{node.id_type.value}->{t_type.value}:{tid}"
                        )
                        queue.append(_Node(tid, t_type, node.path + (t_type.value,)))
                continue

            # At Ensembl gene / gene symbol: fan out to all UniProt (and optional other hops).
            for source in self.sources:
                hits = source.lookup(node.identifier, node.id_type)
                if not hits:
                    continue
                for tid, t_type in hits:
                    attempts.append(
                        f"bridge:{source.name}:{node.id_type.value}->{t_type.value}:{tid}"
                    )
                    if t_type is IdentifierType.UNIPROT or t_type in _BRIDGE_HOP_TYPES:
                        queue.append(_Node(tid, t_type, node.path + (t_type.value,)))

        if not found:
            return [], "", attempts[:50]

        best_path = min(found.values(), key=len)
        return sorted(found.keys()), _format_path(best_path), attempts[:50]


def _match_uniprots_to_entities(
    uniprot_ids: Sequence[str],
    node_type: str,
    graph_entities: set[str],
) -> list[str]:
    """Attach type-specific graph suffixes and keep IDs present in the skeleton."""
    suffix = SUFFIX_NODE_TYPES.get(node_type, "")
    entity_ids: list[str] = []
    for uid in uniprot_ids:
        uid = _strip_isoform(uid)
        candidates = [uid]
        if suffix:
            candidates.append(f"{uid}{suffix}")
        for cand in candidates:
            if cand in graph_entities and cand not in entity_ids:
                entity_ids.append(cand)
    return entity_ids


def _strip_isoform(uniprot_id: str) -> str:
    uid = str(uniprot_id).strip()
    if "-" in uid and uid.split("-")[-1].isdigit():
        return uid.rsplit("-", 1)[0]
    return uid


def _format_path(path: tuple[str, ...]) -> str:
    labels = {
        IdentifierType.GENE_SYMBOL.value: "Gene Symbol",
        IdentifierType.ENSEMBL_GENE.value: "Ensembl Gene",
        IdentifierType.ENSEMBL_TRANSCRIPT.value: "Ensembl Transcript",
        IdentifierType.UNIPROT.value: "UniProt",
    }
    pretty = [labels.get(p, p) for p in path]
    if len(pretty) == 1 and pretty[0] == "UniProt":
        return "Already UniProt"
    return " → ".join(pretty)
