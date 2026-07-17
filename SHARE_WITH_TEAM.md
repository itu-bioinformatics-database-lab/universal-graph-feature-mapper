# Sharing with the team

## What to send

Zip or copy the **entire** `universal_graph_feature_mapper/` directory (~1.9 GB). It already includes:

1. Python package + CLI
2. Universal graph JSON (`data/graphs/`)
3. BioMart / HGNC / UniProt mapper tables (`data/mappers/`)
4. Sample omics for the notebook (`data/sample_omics/`)
5. Notebook tutorial (`notebooks/`)
6. User guide DOCX (`docs/` and Desktop copy)

No separate download of Multi-Omics AD artifacts is required.

## Teammate setup

```bash
cd universal_graph_feature_mapper
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[notebook]"
jupyter notebook notebooks/01_universal_graph_mapping_tutorial.ipynb
```

`config/paths.yaml` already points at packaged `data/` paths.

## Compatibility with Multi-Omics AD

Copy notebook / CLI outputs into the main project:

```
feature_map.csv        → raw/{cohort}/
measured_entities.json → raw/{cohort}/
```

Use the same universal graph version as the rest of the team (this package ships the current `universalGraph_new.json`).
