"""Command-line interface for universal graph feature mapping."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from universal_graph_mapper.config import load_paths
from universal_graph_mapper.constants import ID_COLUMNS
from universal_graph_mapper.biological_graph import BiologicalGraph
from universal_graph_mapper.json_loader import vertex_ids_by_type
from universal_graph_mapper.mapping.pipeline import MappingStrategy, feature_rows_to_frame, run_legacy_then_optimize
from universal_graph_mapper.mapping.prepare import prepare_multi_cohort_raw_with_strategy
from universal_graph_mapper.mapping.resolver import OptimizedResolver
from universal_graph_mapper.mapping.sources import ENHANCEMENT_PROFILES, make_sources_for_profile
from universal_graph_mapper.mapping.stats import compute_dataset_omics_stats, write_mapping_logs, write_unmapped_report


def _parse_omics_files(spec: str) -> dict[str, str]:
    """Parse ``transcript.csv:transcript,protein.csv:protein``."""
    out: dict[str, str] = {}
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(f"Invalid omics spec {part!r}; use filename:node_type")
        fname, ntype = part.split(":", 1)
        out[fname.strip()] = ntype.strip()
    return out


def main_map(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Map omics CSV column names to universal graph node IDs.",
    )
    parser.add_argument("--raw-dir", type=Path, required=True, help="Directory with wide omics CSV files")
    parser.add_argument(
        "--omics",
        required=True,
        help="Comma-separated omics files and graph types, e.g. transcript.csv:transcript,protein.csv:protein",
    )
    parser.add_argument("--graph", type=Path, default=None, help="Universal graph JSON (default: config/paths.yaml)")
    parser.add_argument("--config", type=Path, default=None, help="paths.yaml override")
    parser.add_argument(
        "--strategy",
        choices=["legacy", "optimized", "compare"],
        default="optimized",
        help="legacy = exact match only; optimized = legacy-first then ID DB recovery",
    )
    parser.add_argument(
        "--profile",
        choices=sorted(ENHANCEMENT_PROFILES),
        default="overall",
        help="Offline ID database profile when strategy is optimized",
    )
    parser.add_argument(
        "--override",
        action="append",
        default=[],
        help="Remap omics file to different graph node type, e.g. transcript.csv:gene",
    )
    parser.add_argument("--output-dir", type=Path, default=None, help="Write feature_map.csv here (default: raw-dir)")
    args = parser.parse_args(argv)

    paths = load_paths(args.config)
    graph_path = args.graph or paths.universal_graph
    paths.universal_graph = Path(graph_path)
    paths.validate()

    overrides: dict[str, str] = {}
    for item in args.override:
        fname, ntype = item.split(":", 1)
        overrides[fname.strip()] = ntype.strip()

    out_dir = args.output_dir or args.raw_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    omics_files = _parse_omics_files(args.omics)
    result = prepare_multi_cohort_raw_with_strategy(
        out_dir,
        omics_files,
        graph_path=graph_path,
        omics_read_dir=args.raw_dir if out_dir != args.raw_dir else None,
        omics_type_overrides=overrides or None,
        mapping_strategy=args.strategy,
        mapping_profile=args.profile,
        write_optimized_sidecar=True,
    )

    print(json.dumps({k: v for k, v in result.items() if k != "detailed"}, indent=2, default=str))
    return 0


def main_compare(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compare mapping strategies for one omics file.")
    parser.add_argument("--csv", type=Path, required=True, help="Omics wide CSV")
    parser.add_argument("--node-type", required=True, help="Target graph node type (gene, transcript, protein, ...)")
    parser.add_argument("--graph", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--dataset", default="dataset", help="Label for stats output")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    paths = load_paths(args.config)
    graph_path = args.graph or paths.universal_graph
    paths.universal_graph = Path(graph_path)
    paths.validate()

    graph = BiologicalGraph.from_json_path(graph_path)
    entities = set(vertex_ids_by_type(graph).get(args.node_type, []))
    header = [c for c in pd.read_csv(args.csv, nrows=0).columns if c.strip().lower() not in ID_COLUMNS]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    stats_rows = []
    last_frame = None
    for profile in ("legacy", "biomart", "hgnc", "uniprot", "overall"):
        strategy = MappingStrategy.LEGACY if profile == "legacy" else MappingStrategy.OPTIMIZED
        resolver = (
            OptimizedResolver(sources=make_sources_for_profile(profile), strategy_label=profile)
            if profile != "legacy"
            else OptimizedResolver(sources=[], strategy_label="legacy")
        )
        rows = run_legacy_then_optimize(
            header,
            node_type=args.node_type,
            graph_entities=entities,
            omics_file=args.csv.name,
            resolver=resolver,
            strategy=strategy,
        )
        frame = feature_rows_to_frame(rows)
        last_frame = frame
        frame.to_csv(args.out_dir / f"feature_map_{profile}.csv", index=False)
        stats = compute_dataset_omics_stats(
            frame,
            dataset=args.dataset,
            strategy_label=profile,
            network_nodes_by_omics_file={args.csv.name: len(entities)},
        )
        stats_rows.append(stats)

    stats_df = pd.concat(stats_rows, ignore_index=True)
    stats_df.to_csv(args.out_dir / "coverage_by_strategy.csv", index=False)
    if last_frame is not None:
        write_mapping_logs(last_frame, args.out_dir / args.csv.stem)
        write_unmapped_report(last_frame, args.out_dir / f"unmapped_{args.csv.stem}.tsv")
    print(stats_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "map"
    if cmd == "compare":
        raise SystemExit(main_compare(sys.argv[2:]))
    raise SystemExit(main_map())
