#!/usr/bin/env python3
"""Tests for Ensembl transcript (ENST) resolution via BioMart edge tables."""

from __future__ import annotations

import csv
import json
import unittest
from pathlib import Path

from universal_graph_mapper.mapping.id_detect import IdentifierType, detect_identifier_type
from universal_graph_mapper.mapping.resolver import OptimizedResolver
from universal_graph_mapper.mapping.sources import BioMartSource, _EDGE_TABLES

ROOT = Path(__file__).resolve().parents[1]
BIOMART = ROOT / "data" / "mappers" / "Biomart Mappings"
MAYO_TX = Path(
    "/home/bioinformatics/Desktop/ad-multiomics-graph/benchmark-data/"
    "counts/Mayo_TCX_gene_tx_protein/transcript.csv"
)


class EnsemblTranscriptMappingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = BioMartSource(BIOMART)
        cls.source.load()
        cls.resolver = OptimizedResolver(sources=[cls.source], strategy_label="biomart")

    def test_edge_tables_registered_and_present(self):
        for name in (
            "ensembl_transcript_to_ensembl_gene.json",
            "ensembl_transcript_to_uniprot.json",
            "ensembl_transcript_to_gene_symbol.json",
        ):
            self.assertIn(name, _EDGE_TABLES)
            plain = BIOMART / name
            gz = BIOMART / f"{name}.gz"
            self.assertTrue(plain.is_file() or gz.is_file(), name)
        edges = set(self.source.available_edges())
        self.assertIn((IdentifierType.ENSEMBL_TRANSCRIPT, IdentifierType.ENSEMBL_GENE), edges)
        self.assertIn((IdentifierType.ENSEMBL_TRANSCRIPT, IdentifierType.UNIPROT), edges)
        self.assertIn((IdentifierType.ENSEMBL_TRANSCRIPT, IdentifierType.GENE_SYMBOL), edges)

    def test_detect_enst(self):
        self.assertEqual(
            detect_identifier_type("ENST00000361390"),
            IdentifierType.ENSEMBL_TRANSCRIPT,
        )
        self.assertEqual(
            detect_identifier_type("ENST00000361390.8"),
            IdentifierType.ENSEMBL_TRANSCRIPT,
        )

    def test_resolve_known_coding_transcript(self):
        # MT-ND1 transcript → UniProt P03886 via BioMart
        result = self.resolver.resolve_to_uniprot("ENST00000361390")
        self.assertEqual(result.detected_type, IdentifierType.ENSEMBL_TRANSCRIPT)
        self.assertNotEqual(result.failure_reason, "transcript_without_mapping_table")
        self.assertTrue(result.mapped, result)
        self.assertIn("P03886", result.uniprot_ids)

    def test_mayo_transcript_sample_coverage(self):
        self.assertTrue(MAYO_TX.is_file(), f"missing Mayo transcript CSV: {MAYO_TX}")
        with MAYO_TX.open(newline="") as f:
            cols = next(csv.reader(f))
        feats = [c for c in cols if c not in {"sample_id", "diagnosis"}]
        # Sample systematically across the header for speed.
        sample = feats[:: max(1, len(feats) // 2000)][:2000]
        mapped = 0
        unmapped = 0
        fail_reasons: dict[str, int] = {}
        for fid in sample:
            r = self.resolver.resolve_to_uniprot(fid)
            if r.mapped:
                mapped += 1
            else:
                unmapped += 1
                fail_reasons[r.failure_reason or "unknown"] = (
                    fail_reasons.get(r.failure_reason or "unknown", 0) + 1
                )
        rate = mapped / len(sample)
        # Direct UniProt via transcript or gene hop — expect majority on coding Mayo set.
        self.assertGreaterEqual(
            rate,
            0.35,
            f"mapped={mapped}/{len(sample)} rate={rate:.3f} fails={fail_reasons}",
        )
        # Persist coverage artifact for the commit summary / CI inspection.
        out = ROOT / "outputs" / "ensembl_transcript_mayo_sample_coverage.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            json.dumps(
                {
                    "mayo_file": str(MAYO_TX),
                    "n_sample": len(sample),
                    "n_mapped_to_uniprot": mapped,
                    "n_unmapped": unmapped,
                    "mapping_rate": round(rate, 6),
                    "failure_reasons": fail_reasons,
                },
                indent=2,
            )
            + "\n"
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
