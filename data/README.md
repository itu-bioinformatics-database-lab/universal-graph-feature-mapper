# Packaged data for Universal Graph Feature Mapper

| Path | Contents | In git? |
|------|----------|---------|
| `graphs/universalGraph_new.json` | Universal biological graph (~1.3 GB) | **No** (gitignored) |
| `mappers/{Biomart,HGNC,UniProt} Mappings/` | Offline ID edge tables | Yes (except `master_identifier_mappings_long.json` >100 MB) |
| `sample_omics/` | Tutorial CSVs | Yes |

## Obtaining the universal graph

```bash
cp /path/to/Multi-Omics-AD/artifacts/graphs/ROSMAP/universalGraph_new.json \
   data/graphs/universalGraph_new.json
```

`config/paths.yaml` already points at `data/graphs/universalGraph_new.json`.
