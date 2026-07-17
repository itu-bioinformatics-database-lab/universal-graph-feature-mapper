# Universal Graph Feature Mapper

Standalone, **self-contained** package for mapping omics feature columns onto a universal biological graph (MOGI-compatible JSON).

All runtime dependencies are under `data/` — no Multi-Omics AD checkout required.

## Contents

| Path | Role |
|------|------|
| `data/graphs/universalGraph_new.json` | Universal graph (~1.3 GB) |
| `data/mappers/{Biomart,HGNC,UniProt} Mappings/` | Offline ID tables (~600 MB) |
| `data/sample_omics/` | Tutorial CSV matrices |
| `notebooks/01_universal_graph_mapping_tutorial.ipynb` | Full usage tutorial |
| `config/paths.yaml` | Default paths (relative to this package) |
| `docs/Universal_Graph_Feature_Mapper_User_Guide.docx` | Detailed user guide |

## Install

```bash
cd universal_graph_feature_mapper
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[notebook]"
```

## Notebook tutorial

```bash
jupyter notebook notebooks/01_universal_graph_mapping_tutorial.ipynb
```

The notebook validates packaged data, maps mixed IDs, prepares transcript+protein feature maps, and compares legacy vs BioMart/HGNC/UniProt/overall strategies. Outputs go to `outputs/notebook_tutorial/`.

## CLI

```bash
ugm-map \
  --raw-dir data/sample_omics \
  --omics "transcript.csv:transcript,protein.csv:protein" \
  --override "transcript.csv:gene" \
  --strategy optimized \
  --profile overall \
  --output-dir outputs/cli_prepared
```

## Python API

```python
from pathlib import Path
from universal_graph_mapper.mapping.prepare import prepare_multi_cohort_raw_with_strategy

prepare_multi_cohort_raw_with_strategy(
    Path("outputs/my_cohort"),
    omics_files={"transcript.csv": "transcript", "protein.csv": "protein"},
    omics_read_dir=Path("data/sample_omics"),
    omics_type_overrides={"transcript.csv": "gene"},
    mapping_strategy="optimized",
    mapping_profile="overall",
)
```

## Graph node ID conventions

| Node type | Vertex ID | Example |
|-----------|-----------|---------|
| gene | bare UniProt | `Q07973` |
| transcript | accession + `_transcript` | `Q07973_transcript` |
| protein | accession + `_protein` | `Q07973_protein` |

## Sharing

Zip or rsync this entire folder (~1.9 GB). See `SHARE_WITH_TEAM.md`.
