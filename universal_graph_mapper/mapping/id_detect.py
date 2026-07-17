"""Automatic biological identifier type detection."""

from __future__ import annotations

import re
from enum import Enum


class IdentifierType(str, Enum):
    UNIPROT = "uniprot"
    ENSEMBL_GENE = "ensembl_gene"
    ENSEMBL_TRANSCRIPT = "ensembl_transcript"
    GENE_SYMBOL = "gene_symbol"
    OTHER = "other"
    EMPTY = "empty"


# UniProtKB accession (including A0A* isoform bases).
_UNIPROT_RE = re.compile(
    r"^(?:[A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[OPQ][0-9][A-Z0-9]{3}[0-9]|A0A[A-Z0-9]{7}"
    r")(?:-\d+)?$",
    re.IGNORECASE,
)
_ENSEMBL_GENE_RE = re.compile(r"^ENSG\d{11}(?:\.\d+)?$", re.IGNORECASE)
_ENSEMBL_TX_RE = re.compile(r"^ENST\d{11}(?:\.\d+)?$", re.IGNORECASE)
_GRAPH_SUFFIX_RE = re.compile(r"_(?:transcript|protein)$", re.IGNORECASE)
_VERSION_RE = re.compile(r"\.\d+$")


def strip_graph_suffix(identifier: str) -> str:
    """Remove project graph suffixes ``_transcript`` / ``_protein``."""
    return _GRAPH_SUFFIX_RE.sub("", str(identifier).strip())


def strip_version(identifier: str) -> str:
    """Strip Ensembl version suffix (``.12``)."""
    return _VERSION_RE.sub("", identifier)


def normalize_identifier(identifier: str) -> str:
    """Normalize for BioMart lookup (strip graph suffix + Ensembl version)."""
    raw = strip_graph_suffix(str(identifier).strip())
    return strip_version(raw)


def detect_identifier_type(identifier: str) -> IdentifierType:
    """Infer identifier type for mixed omics columns.

    Never assumes a single type for a whole file; call per feature.
    """
    if identifier is None:
        return IdentifierType.EMPTY
    raw = str(identifier).strip()
    if not raw:
        return IdentifierType.EMPTY

    base = normalize_identifier(raw)
    if not base:
        return IdentifierType.EMPTY
    if _UNIPROT_RE.match(base):
        return IdentifierType.UNIPROT
    if _ENSEMBL_TX_RE.match(base):
        return IdentifierType.ENSEMBL_TRANSCRIPT
    if _ENSEMBL_GENE_RE.match(base):
        return IdentifierType.ENSEMBL_GENE
    # Gene symbols: alphanumeric with limited punctuation (TP53, MT-ND1, HLA-A).
    if re.match(r"^[A-Za-z][A-Za-z0-9._@/-]*$", base) and not base.startswith("ENS"):
        return IdentifierType.GENE_SYMBOL
    return IdentifierType.OTHER
