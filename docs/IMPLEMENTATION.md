# Implementation Strategy — Universal Graph Feature Mapper

This document explains **how** the package maps omics feature identifiers onto a universal
biological graph. For usage instructions see `README.md` and the DOCX user guide.

---

## 1. Problem

Given omics matrices whose **columns** are arbitrary feature identifiers
(Ensembl gene IDs, HGNC gene symbols, UniProt accessions, or already-suffixed graph IDs),
attach each column to a **vertex ID in the universal graph** so measurements can be placed
on the correct node.

The universal graph is UniProt-centric:

| Node type | Vertex ID form | Example |
|-----------|----------------|---------|
| gene | bare UniProt accession | `Q07973` |
| transcript | accession + `_transcript` | `Q07973_transcript` |
| protein | accession + `_protein` | `Q07973_protein` |

Suffix rules live in `graph_schema.SUFFIX_NODE_TYPES`.

---

## 2. Two-stage strategy (legacy-first, then optimized recovery)

Mapping is intentionally **layered**. The controlling invariant is:

> **A valid legacy (exact) hit is never overwritten by optimized recovery.**

### Stage A — Legacy exact match (`pipeline.legacy_map_column`)

For each column, candidates are tried in order:

1. the raw column name,
2. the stripped column name,
3. the name with the type suffix applied (`apply_entity_suffix`),
4. for `gene` targets, a `_transcript`-stripped base (transcript→gene override case).

The first candidate that exists in the graph's entity set wins. This is fast, deterministic,
and covers columns that already name a graph vertex.

### Stage B — Optimized recovery (`resolver.OptimizedResolver`)

Applied **only** to columns that legacy could not map, and **only** for the optimizable node
types `{gene, transcript, protein}` (`OPTIMIZABLE_NODE_TYPES`). Everything else is passed
through as a legacy result.

Strategy selection (`MappingStrategy`):

- `legacy` — Stage A only.
- `optimized` — Stage A, then Stage B on the remainder (**default**).
- `compare` — behaves like `optimized` for the written map, but also emits per-strategy
  frames for coverage comparison.

---

## 3. Identifier detection (`id_detect.py`)

Each column is typed **independently** (a file is never assumed to be uniform):

1. empty → `EMPTY`
2. UniProtKB regex (incl. `A0A…`) → `UNIPROT`
3. `ENST…` → `ENSEMBL_TRANSCRIPT` (currently unrecoverable — no ENST edge tables)
4. `ENSG…` → `ENSEMBL_GENE`
5. symbol-like token → `GENE_SYMBOL`
6. otherwise → `OTHER`

Normalization (`normalize_identifier`) strips graph suffixes (`_transcript`, `_protein`)
and Ensembl version suffixes (`.12`); isoforms like `P04637-2` reduce to `P04637`.

---

## 4. Offline mapping sources (`sources.py`)

Recovery uses **offline JSON edge tables** — no network calls. Each source is a set of
one-hop lookup tables sharing a common record schema:

```json
{"source_id": "TP53", "source_id_type": "gene_symbol",
 "target_id": "P04637", "target_id_type": "uniprot_id",
 "mapping_source": "ensembl_biomart"}
```

Three interchangeable sources, each a directory of tables:

| Source | Directory (config) |
|--------|--------------------|
| `BioMartSource` | `data/mappers/Biomart Mappings/` |
| `HgncSource` | `data/mappers/HGNC Mappings/` |
| `UniProtIdMappingSource` | `data/mappers/UniProt Mappings/` |

Profiles (`ENHANCEMENT_PROFILES`): `biomart`, `hgnc`, `uniprot`, and `overall`
(all three; the BFS explores every source at each hop). `overall` is the recommended default.

One-to-many relationships (one symbol/Ensembl → many UniProts) are **preserved**, not collapsed.

---

## 5. Multi-hop resolution (`resolver.resolve_to_uniprot`)

For a non-UniProt identifier, the resolver runs a **breadth-first search** over the plugged-in
sources until it reaches UniProt accessions:

```
feature id  →  (BioMart/HGNC/UniProt hops)  →  {UniProt accessions}
```

- Bounded by `max_depth` (default 4).
- Visited `(id, type)` pairs are de-duplicated.
- Already-UniProt inputs short-circuit (isoform stripped).
- `ENSEMBL_TRANSCRIPT` and unrecognized types fail fast with an explicit `failure_reason`.
- The discovery path is recorded for provenance (`mapping_path`, e.g. `Ensembl Gene → UniProt`).

---

## 6. Graph attachment + UniProt sibling bridge (`resolve_to_graph_entities`)

Resolved UniProt accessions are attached to graph vertices via
`_match_uniprots_to_entities` (applies the node-type suffix and keeps only IDs present in
the graph).

If **none** of the direct UniProts are in the graph, the resolver attempts a **sibling
bridge**:

```
UniProtB (in data, NOT in graph)
   → Ensembl gene / gene symbol         (bridge hop types)
   → {UniProtA, UniProtB, ...}          (fan-out)
   → keep UniProtA if it IS in the graph
```

This recovers synonymous in-network accessions instead of dropping a feature just because the
data's specific accession is absent from this graph build. Bridge hops are limited to
`ENSEMBL_GENE` and `GENE_SYMBOL` (`_BRIDGE_HOP_TYPES`). Original out-of-network IDs are kept
in the provenance list; the mapping is flagged `uniprot_not_in_network` only if the bridge
also fails.

---

## 7. Orchestration (`pipeline.run_legacy_then_optimize`)

Per column the pipeline emits a `FeatureMapRow` with full provenance:

| Field | Meaning |
|-------|---------|
| `omics_file`, `feature_col` | source location |
| `entity_id`, `entity_ids` | primary + all matched graph vertex IDs |
| `node_type` | target graph type (after any override) |
| `in_graph` | whether a graph vertex was matched |
| `strategy` | `legacy` or the recovery profile that succeeded |
| `detected_type` | inferred identifier type |
| `mapping_path` | human-readable hop path |
| `uniprot_ids` | intermediate UniProt accessions |
| `status`, `failure_reason` | success/unmapped + diagnostic |

---

## 8. Dataset preparation (`prepare.prepare_multi_cohort_raw_with_strategy`)

The top-level entry point:

1. Loads the universal graph (from `config/paths.yaml` if `graph_path` is omitted).
2. For each omics file, resolves the target node type (honoring `omics_type_overrides`,
   e.g. `transcript.csv → gene`).
3. Reads only the header (feature columns) from each matrix; ID columns
   (`sample_id`, `diagnosis`, `label`, `label_id`) are ignored.
4. Runs legacy-then-optimized mapping per column.
5. Writes outputs:
   - `feature_map.csv` — canonical `omics_file, feature_col, entity_id, node_type, in_graph`
   - `feature_map_detailed.csv` — full provenance sidecar
   - `measured_entities.json` — sorted unique matched graph entity IDs (seeds for downstream
     graph pruning / model building)

`omics_read_dir` lets you read matrices from one folder while writing maps to another.

---

## 9. Coverage rates (`stats.py`)

Two rates quantify success:

- **From data** = mapped features ÷ input features (recall over the dataset).
- **To network** = unique graph vertices hit ÷ graph vertices of that type (coverage of the
  network).

The `ugm-compare` CLI and notebook Section 9 report both, per strategy, so legacy vs optimized
gains are explicit.

---

## 10. Design invariants (summary)

1. **Legacy-first** — exact graph hits are authoritative and never overwritten.
2. **No invented nodes** — mapping only ever attaches to vertices that already exist.
3. **One-to-many preserved** — a single symbol/Ensembl may map to multiple graph vertices.
4. **Per-column typing** — identifier type is detected per feature, not per file.
5. **Suffix contract** — `_transcript` / `_protein` suffixes are part of the graph ID scheme.
6. **Offline & deterministic** — all recovery uses shipped JSON tables; no network calls.
7. **Provenance kept** — every mapping records how it was derived and why it failed.

---

## 11. Known limitations

- **ENST (Ensembl transcript) IDs** are not recoverable — no ENST edge tables are shipped.
- **Metabolites** have no vertex type in the universal graph and are out of scope.
- Coverage depends on the **specific graph build**; a different `universalGraph_new.json`
  version can change which accessions are in-network.

---

## 12. Code map

| Module | Responsibility |
|--------|----------------|
| `mapping/id_detect.py` | identifier typing + normalization |
| `mapping/sources.py` | offline BioMart / HGNC / UniProt edge tables + profiles |
| `mapping/resolver.py` | BFS to UniProt, graph match, sibling bridge |
| `mapping/pipeline.py` | legacy-first orchestration, `FeatureMapRow` |
| `mapping/prepare.py` | dataset-level prepare + output files |
| `mapping/stats.py` | coverage tables / rates |
| `graph_schema.py` | suffix rules, CSV↔type mapping |
| `biological_graph.py`, `json_loader.py` | graph loading + vertex grouping |
| `config.py` | path resolution (`paths.yaml` + `UGM_*` env) |
| `cli.py` | `ugm-map`, `ugm-compare` |
