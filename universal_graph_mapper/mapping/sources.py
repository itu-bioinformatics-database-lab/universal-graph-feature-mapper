"""Pluggable offline mapping sources (BioMart, HGNC, UniProt ID mapping)."""

from __future__ import annotations

import gzip
import json
from abc import ABC, abstractmethod
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Iterable, Mapping, Sequence

from universal_graph_mapper.mapping.id_detect import IdentifierType

def _default_dirs():
    from universal_graph_mapper.config import load_paths
    p = load_paths()
    return p.biomart_dir, p.hgnc_dir, p.uniprot_dir


def DEFAULT_BIOMART_DIR() -> Path:
    return _default_dirs()[0]


def DEFAULT_HGNC_DIR() -> Path:
    return _default_dirs()[1]


def DEFAULT_UNIPROT_DIR() -> Path:
    return _default_dirs()[2]

# Shared list-of-record JSON edge tables across BioMart / HGNC / UniProt folders.
_EDGE_TABLES = (
    "gene_symbol_to_uniprot.json",
    "ensembl_gene_to_uniprot.json",
    "gene_symbol_to_ensembl_gene.json",
    "ensembl_gene_to_gene_symbol.json",
    "uniprot_to_ensembl_gene.json",
    "uniprot_to_gene_symbol.json",
    # Ensembl transcript edges (BioMart; absent folders simply skip missing files)
    "ensembl_transcript_to_ensembl_gene.json",
    "ensembl_transcript_to_uniprot.json",
    "ensembl_transcript_to_gene_symbol.json",
)


class MappingSource(ABC):
    """One hop: source ID → one or more target IDs."""

    name: str = "abstract"

    @abstractmethod
    def lookup(self, source_id: str, source_type: IdentifierType) -> list[tuple[str, IdentifierType]]:
        """Return (target_id, target_type) pairs. Empty if unknown."""


class OfflineJsonMappingSource(MappingSource):
    """Offline JSON edge tables (one-to-many preserved). Shared schema for all DBs."""

    def __init__(
        self,
        tables_dir: str | Path,
        *,
        name: str,
        table_names: Sequence[str] = _EDGE_TABLES,
    ):
        self.tables_dir = Path(tables_dir)
        self.name = name
        self.table_names = tuple(table_names)
        # (source_type, target_type) → source_id → [target_ids]
        self._maps: dict[tuple[IdentifierType, IdentifierType], DefaultDict[str, list[str]]] = {}
        self._loaded = False

    def load(self) -> None:
        if self._loaded:
            return
        if not self.tables_dir.is_dir():
            raise FileNotFoundError(f"Mapping directory not found: {self.tables_dir}")
        for fname in self.table_names:
            path = self._resolve_table_path(fname)
            if path is None:
                continue
            records = _load_json_records(path)
            if not isinstance(records, list) or not records:
                continue
            s_type = _id_type(records[0].get("source_id_type", ""))
            t_type = _id_type(records[0].get("target_id_type", ""))
            if s_type is None or t_type is None:
                continue
            bucket = self._maps.setdefault((s_type, t_type), defaultdict(list))
            for row in records:
                sid = str(row.get("source_id", "")).strip()
                tid = str(row.get("target_id", "")).strip()
                if not sid or not tid:
                    continue
                if tid not in bucket[sid]:
                    bucket[sid].append(tid)
        self._loaded = True

    def _resolve_table_path(self, fname: str) -> Path | None:
        """Prefer plain JSON; fall back to ``.json.gz`` for large BioMart tables."""
        path = self.tables_dir / fname
        if path.exists():
            return path
        if fname.endswith(".json"):
            gz = self.tables_dir / f"{fname}.gz"
            if gz.exists():
                return gz
        return None

    def lookup(self, source_id: str, source_type: IdentifierType) -> list[tuple[str, IdentifierType]]:
        self.load()
        out: list[tuple[str, IdentifierType]] = []
        for (s_type, t_type), table in self._maps.items():
            if s_type != source_type:
                continue
            for tid in table.get(source_id, []):
                out.append((tid, t_type))
        return out

    def available_edges(self) -> list[tuple[IdentifierType, IdentifierType]]:
        self.load()
        return sorted(self._maps.keys(), key=lambda x: (x[0].value, x[1].value))


class BioMartSource(OfflineJsonMappingSource):
    """Offline Ensembl BioMart JSON tables."""

    def __init__(self, tables_dir: str | Path | None = None):
        super().__init__(tables_dir or DEFAULT_BIOMART_DIR(), name="biomart")


class HgncSource(OfflineJsonMappingSource):
    """Offline HGNC JSON tables."""

    def __init__(self, tables_dir: str | Path | None = None):
        super().__init__(tables_dir or DEFAULT_HGNC_DIR(), name="hgnc")


class UniProtIdMappingSource(OfflineJsonMappingSource):
    """Offline UniProt ID-mapping JSON tables."""

    def __init__(self, tables_dir: str | Path | None = None):
        super().__init__(tables_dir or DEFAULT_UNIPROT_DIR(), name="uniprot")


class DictMappingSource(MappingSource):
    """In-memory / test source: mapping[(source_type, target_type)][id] = [targets]."""

    name = "dict"

    def __init__(
        self,
        tables: Mapping[tuple[IdentifierType, IdentifierType], Mapping[str, Iterable[str]]],
        *,
        name: str = "dict",
    ):
        self.name = name
        self._maps = {
            edge: {k: list(v) for k, v in table.items()} for edge, table in tables.items()
        }

    def lookup(self, source_id: str, source_type: IdentifierType) -> list[tuple[str, IdentifierType]]:
        out: list[tuple[str, IdentifierType]] = []
        for (s_type, t_type), table in self._maps.items():
            if s_type != source_type:
                continue
            for tid in table.get(source_id, []):
                out.append((str(tid), t_type))
        return out


# Enhancement profiles used by comparison / prepare (legacy is applied first separately).
ENHANCEMENT_PROFILES: dict[str, tuple[str, ...]] = {
    "biomart": ("biomart",),
    "hgnc": ("hgnc",),
    "uniprot": ("uniprot",),
    # Overall: try all DBs; BFS explores every source at each hop.
    "overall": ("biomart", "hgnc", "uniprot"),
}


def make_source(name: str, tables_dir: str | Path | None = None) -> MappingSource:
    """Build a named offline mapping source."""
    key = str(name).strip().lower()
    if key in {"biomart", "ensembl_biomart", "ensembl"}:
        return BioMartSource(tables_dir)
    if key in {"hgnc"}:
        return HgncSource(tables_dir)
    if key in {"uniprot", "uniprot_idmapping", "uniprot_mapping"}:
        return UniProtIdMappingSource(tables_dir)
    raise ValueError(f"Unknown mapping source: {name!r}")


def make_sources_for_profile(profile: str) -> list[MappingSource]:
    """Return MappingSource instances for an enhancement profile key."""
    key = str(profile).strip().lower()
    names = ENHANCEMENT_PROFILES.get(key)
    if names is None:
        raise ValueError(
            f"Unknown enhancement profile {profile!r}; "
            f"expected one of {sorted(ENHANCEMENT_PROFILES)}"
        )
    return [make_source(n) for n in names]


def _id_type(label: str) -> IdentifierType | None:
    key = str(label).strip().lower()
    return {
        "uniprot_id": IdentifierType.UNIPROT,
        "uniprot": IdentifierType.UNIPROT,
        "ensembl_gene_id": IdentifierType.ENSEMBL_GENE,
        "ensembl_gene": IdentifierType.ENSEMBL_GENE,
        "ensembl_transcript_id": IdentifierType.ENSEMBL_TRANSCRIPT,
        "ensembl_transcript": IdentifierType.ENSEMBL_TRANSCRIPT,
        "gene_symbol": IdentifierType.GENE_SYMBOL,
        "hgnc_symbol": IdentifierType.GENE_SYMBOL,
    }.get(key)


def _load_json_records(path: Path):
    """Load a JSON list from ``.json`` or ``.json.gz``."""
    if path.suffix == ".gz" or str(path).endswith(".json.gz"):
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            return json.load(fh)
    return json.loads(path.read_text(encoding="utf-8"))
