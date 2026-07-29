"""Mapping statistics, comparison tables, logs, and report generation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def compute_dataset_omics_stats(
    mapped_df: pd.DataFrame,
    *,
    dataset: str,
    strategy_label: str,
    network_nodes_by_omics_file: dict[str, int],
) -> pd.DataFrame:
    """Per-omics counts and dual rates (data-side and network-side).

    Rates
    -----
    mapping_rate_from_data = |intersection| / |data features|
        (features with ≥1 in-graph entity) / (total input features)

    mapping_rate_to_network = |intersection nodes| / |graph nodes|
        (unique graph entities hit) / (nodes of that type in the network)
    """
    rows = []
    if mapped_df.empty:
        return pd.DataFrame()
    for omics_file, grp in mapped_df.groupby("omics_file"):
        n_data = len(grp)
        mapped_grp = grp[grp["in_graph"].astype(bool)]
        n_mapped_features = int(len(mapped_grp))
        unique_entities: set[str] = set()
        if "entity_ids" in mapped_grp.columns:
            for cell in mapped_grp["entity_ids"].fillna("").astype(str):
                for eid in cell.split("|"):
                    eid = eid.strip()
                    if eid:
                        unique_entities.add(eid)
        if "entity_id" in mapped_grp.columns:
            for eid in mapped_grp["entity_id"].fillna("").astype(str):
                eid = eid.strip()
                if eid:
                    unique_entities.add(eid)
        n_intersection_nodes = len(unique_entities)
        n_graph = int(network_nodes_by_omics_file.get(str(omics_file), 0))
        node_type = grp["node_type"].iloc[0] if len(grp) else ""
        rate_from_data = n_mapped_features / max(n_data, 1)
        rate_to_network = n_intersection_nodes / max(n_graph, 1)

        rows.append(
            {
                "dataset": dataset,
                "omics_file": omics_file,
                "omics": _omics_label(str(omics_file), node_type),
                "node_type": node_type,
                "strategy": strategy_label,
                "n_graph_nodes": n_graph,
                "n_data_features": n_data,
                "n_mapped_features": n_mapped_features,
                "n_intersection_nodes": n_intersection_nodes,
                "n_unmapped_features": n_data - n_mapped_features,
                "n_unique_unmapped_features": int(
                    grp.loc[~grp["in_graph"].astype(bool), "feature_col"].nunique()
                ),
                # Aliases kept for older report code.
                "n_input_features": n_data,
                "feature_mapping_rate": rate_from_data,
                "mapping_rate_from_data": rate_from_data,
                "mapping_rate_to_network": rate_to_network,
                "n_unique_uniprot_nodes": len({_base_uniprot(e) for e in unique_entities}),
                "n_unique_graph_entities": n_intersection_nodes,
                "network_protein_nodes": n_graph,
                "network_coverage": rate_to_network,
            }
        )
    return pd.DataFrame(rows)


def comparison_table(legacy_stats: pd.DataFrame, optimized_stats: pd.DataFrame) -> pd.DataFrame:
    """Side-by-side legacy vs one enhancement (back-compat; prefixes legacy_/optimized_)."""
    return multi_strategy_comparison_table(
        {"legacy": legacy_stats, "optimized": optimized_stats},
        baseline="legacy",
        primary="optimized",
    )


def multi_strategy_comparison_table(
    stats_by_strategy: dict[str, pd.DataFrame],
    *,
    baseline: str = "legacy",
    primary: str | None = None,
) -> pd.DataFrame:
    """Wide table: absolute counts + rates for every strategy vs a baseline.

    Column naming (for strategy ``S``)::

        {S}_mapped_features, {S}_intersection_nodes,
        {S}_rate_from_data, {S}_rate_to_network, ...

    Deltas are baseline → each non-baseline strategy (and an extra optimized_*
    alias when ``primary`` is set for plot back-compat).
    """
    if not stats_by_strategy:
        return pd.DataFrame()
    keys = ["dataset", "omics_file", "omics", "node_type"]
    size_cols = ["n_graph_nodes", "n_data_features"]

    frames = [(k, df) for k, df in stats_by_strategy.items() if df is not None and not df.empty]
    if not frames:
        return pd.DataFrame()

    # Union of size keys from all strategies.
    left = frames[0][1][keys + size_cols].copy()
    for _, df in frames[1:]:
        left = keys_merge_sizes(left, df, keys, size_cols)

    merged = left.copy()
    for label, df in frames:
        prefix = f"{label}_"
        merged = merged.merge(_prefix_strategy_stats(df, keys, prefix), on=keys, how="outer")

    base = baseline if baseline in stats_by_strategy else frames[0][0]
    # Back-compat aliases when comparing exactly legacy vs optimized.
    if "legacy" in stats_by_strategy and "optimized" in stats_by_strategy:
        merged["legacy_feature_mapping_rate"] = merged.get("legacy_rate_from_data")
        merged["optimized_feature_mapping_rate"] = merged.get("optimized_rate_from_data")
        merged["legacy_network_coverage"] = merged.get("legacy_rate_to_network")
        merged["optimized_network_coverage"] = merged.get("optimized_rate_to_network")

    for label, _ in frames:
        if label == base:
            continue
        for metric, out_name in (
            ("mapped_features", f"delta_{label}_mapped_features"),
            ("intersection_nodes", f"delta_{label}_intersection_nodes"),
            ("unique_uniprot", f"delta_{label}_unique_uniprot"),
            ("rate_from_data", f"delta_{label}_rate_from_data"),
            ("rate_to_network", f"delta_{label}_rate_to_network"),
        ):
            bcol = f"{base}_{metric}"
            scol = f"{label}_{metric}"
            if bcol in merged.columns and scol in merged.columns:
                merged[out_name] = merged[scol].fillna(0) - merged[bcol].fillna(0)

    # Primary enhancement deltas (optimized_* / delta_* without strategy infix).
    prim = primary
    if prim is None:
        for candidate in ("overall", "optimized", "biomart"):
            if candidate in stats_by_strategy and candidate != base:
                prim = candidate
                break
    if prim and prim in stats_by_strategy and prim != base:
        for metric in (
            "mapped_features",
            "intersection_nodes",
            "unique_uniprot",
            "rate_from_data",
            "rate_to_network",
        ):
            alt = {
                "mapped_features": "delta_mapped_features",
                "intersection_nodes": "delta_intersection_nodes",
                "unique_uniprot": "delta_unique_uniprot",
                "rate_from_data": "delta_feature_mapping_rate",
                "rate_to_network": "delta_rate_to_network",
            }[metric]
            src = f"delta_{prim}_{metric}"
            if src in merged.columns:
                merged[alt] = merged[src]
        merged["delta_rate_from_data"] = merged.get("delta_feature_mapping_rate")
        merged["delta_network_coverage"] = merged.get("delta_rate_to_network")
        # Alias primary as optimized_* for existing bar plots when primary != optimized.
        if prim != "optimized":
            for metric in (
                "mapped_features",
                "intersection_nodes",
                "unique_uniprot",
                "rate_from_data",
                "rate_to_network",
            ):
                src = f"{prim}_{metric}"
                if src in merged.columns:
                    merged[f"optimized_{metric}"] = merged[src]
            merged["optimized_feature_mapping_rate"] = merged.get(f"{prim}_rate_from_data")
            merged["optimized_network_coverage"] = merged.get(f"{prim}_rate_to_network")

    return merged


def _prefix_strategy_stats(df: pd.DataFrame, keys: list[str], prefix: str) -> pd.DataFrame:
    keep = keys + [
        "n_mapped_features",
        "n_intersection_nodes",
        "mapping_rate_from_data",
        "mapping_rate_to_network",
        "n_unique_uniprot_nodes",
        "n_unique_graph_entities",
    ]
    keep = [c for c in keep if c in df.columns]
    out = df[keep].copy()
    out = out.rename(columns={c: f"{prefix}{c}" for c in keep if c not in keys})
    out = out.rename(
        columns={
            f"{prefix}n_mapped_features": f"{prefix}mapped_features",
            f"{prefix}n_intersection_nodes": f"{prefix}intersection_nodes",
            f"{prefix}mapping_rate_from_data": f"{prefix}rate_from_data",
            f"{prefix}mapping_rate_to_network": f"{prefix}rate_to_network",
            f"{prefix}n_unique_uniprot_nodes": f"{prefix}unique_uniprot",
            f"{prefix}n_unique_graph_entities": f"{prefix}unique_entities",
        }
    )
    return out


def format_coverage_table(
    comparison: pd.DataFrame,
    *,
    strategies: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Human-facing absolute counts + both rates for legacy and enhancements."""
    if comparison.empty:
        return comparison
    if strategies is None:
        # Infer strategy prefixes from column names like `{s}_mapped_features`.
        found: list[str] = []
        for c in comparison.columns:
            if c.endswith("_mapped_features") and not c.startswith("delta_"):
                found.append(c[: -len("_mapped_features")])
        # Prefer a stable order.
        preferred = ["legacy", "biomart", "hgnc", "uniprot", "overall", "optimized"]
        strategies = [s for s in preferred if s in found] + [s for s in found if s not in preferred]

    cols = ["dataset", "omics", "node_type", "n_graph_nodes", "n_data_features"]
    for s in strategies:
        cols.extend(
            [
                f"{s}_mapped_features",
                f"{s}_intersection_nodes",
                f"{s}_rate_from_data",
                f"{s}_rate_to_network",
            ]
        )
    # Prefer overall (or optimized) deltas when present.
    for delta_col in (
        "delta_overall_mapped_features",
        "delta_mapped_features",
        "delta_overall_intersection_nodes",
        "delta_intersection_nodes",
    ):
        if delta_col in comparison.columns and delta_col not in cols:
            cols.append(delta_col)

    available = [c for c in cols if c in comparison.columns]
    out = comparison[available].copy()
    for c in available:
        if "rate" in c:
            out[c] = out[c].map(lambda v: f"{100 * float(v):.2f}%" if pd.notna(v) else "")
        elif c.startswith("n_") or "mapped" in c or "intersection" in c or c.startswith("delta_"):
            out[c] = out[c].map(lambda v: int(v) if pd.notna(v) else 0)
    return out


def plot_strategies_grouped_bars(
    comparison: pd.DataFrame,
    *,
    metric: str,
    strategies: Sequence[str],
    title: str,
    ylabel: str,
    out_path: Path,
    labels: Sequence[str] | None = None,
) -> None:
    """Grouped bars for one metric across multiple mapping strategies."""
    if comparison.empty:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    row_labels = labels or [f"{r.dataset}\n{r.omics}" for r in comparison.itertuples()]
    x = np.arange(len(row_labels))
    n = len(strategies)
    width = min(0.8 / max(n, 1), 0.2)
    fig, ax = plt.subplots(figsize=(max(9, len(row_labels) * 1.3), 5))
    display = {
        "legacy": "Legacy",
        "biomart": "BioMart",
        "hgnc": "HGNC",
        "uniprot": "UniProt",
        "overall": "Overall",
        "optimized": "Optimized",
    }
    for i, strat in enumerate(strategies):
        col = f"{strat}_{metric}"
        if col not in comparison.columns:
            continue
        offset = (i - (n - 1) / 2) * width
        ax.bar(x + offset, comparison[col].fillna(0), width, label=display.get(strat, strat))
    ax.set_xticks(x)
    ax.set_xticklabels(row_labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)



def keys_merge_sizes(
    legacy_stats: pd.DataFrame,
    optimized_stats: pd.DataFrame,
    keys: list[str],
    size_cols: list[str],
) -> pd.DataFrame:
    """Take size columns from either frame (prefer max / non-null)."""
    a = legacy_stats[keys + size_cols].copy()
    b = optimized_stats[keys + size_cols].copy()
    merged = a.merge(b, on=keys, how="outer", suffixes=("_l", "_r"))
    out = merged[keys].copy()
    for c in size_cols:
        out[c] = merged[[f"{c}_l", f"{c}_r"]].bfill(axis=1).iloc[:, 0]
    return out


def write_mapping_logs(mapped_df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "feature_col",
        "detected_type",
        "strategy",
        "mapping_path",
        "uniprot_ids",
        "entity_id",
        "status",
        "failure_reason",
    ]
    use = [c for c in cols if c in mapped_df.columns]
    mapped_df[use].to_csv(out_path, sep="\t", index=False)


def write_unmapped_report(mapped_df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    unmapped = mapped_df[~mapped_df["in_graph"].astype(bool)].copy()
    unmapped.to_csv(out_path, sep="\t", index=False)


def plot_legacy_vs_optimized_bars(
    comparison: pd.DataFrame,
    value_col_legacy: str,
    value_col_opt: str,
    *,
    title: str,
    ylabel: str,
    out_path: Path,
) -> None:
    if comparison.empty:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    labels = [
        f"{r.dataset}\n{r.omics}" for r in comparison.itertuples()
    ]
    x = np.arange(len(labels))
    width = 0.35
    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 1.2), 5))
    ax.bar(x - width / 2, comparison[value_col_legacy].fillna(0), width, label="Legacy")
    ax.bar(x + width / 2, comparison[value_col_opt].fillna(0), width, label="Optimized")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_mapping_sources_pie(mapped_df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    paths = mapped_df["mapping_path"].fillna("Not mapped").astype(str)
    counts = paths.value_counts()
    # Collapse long tails for readability.
    if len(counts) > 8:
        top = counts.iloc[:7]
        other = counts.iloc[7:].sum()
        counts = pd.concat([top, pd.Series({"Other": other})])
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.pie(counts.values, labels=counts.index, autopct="%1.1f%%", startangle=90)
    ax.set_title("Mapping Sources / Paths")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_recovery(
    comparison: pd.DataFrame,
    out_path: Path,
) -> None:
    """Original → recovered → final mapped feature counts."""
    if comparison.empty:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(max(8, len(comparison) * 1.1), 5))
    labels = [f"{r.dataset}/{r.omics}" for r in comparison.itertuples()]
    legacy = comparison["legacy_mapped_features"].fillna(0).to_numpy()
    delta = comparison["delta_mapped_features"].clip(lower=0).fillna(0).to_numpy()
    final = comparison["optimized_mapped_features"].fillna(0).to_numpy()
    x = np.arange(len(labels))
    ax.bar(x, legacy, label="Original (legacy)")
    ax.bar(x, delta, bottom=legacy, label="Recovered (optimized)")
    ax.plot(x, final, "ko-", label="Final mapped")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Mapped features")
    ax.set_title("Mapping Recovery")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_sankey_flows(mapped_df: pd.DataFrame, out_path: Path) -> None:
    """Approximate Sankey via stacked horizontal flow counts (no external sankey dep)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if mapped_df.empty or "detected_type" not in mapped_df.columns:
        return
    stages = [
        ("detected_type", "Detected type"),
        ("mapping_path", "Mapping path"),
        ("status", "Status"),
    ]
    fig, axes = plt.subplots(1, len(stages), figsize=(14, 5))
    for ax, (col, title) in zip(axes, stages):
        counts = mapped_df[col].fillna("unknown").astype(str).value_counts().head(10)
        ax.barh(counts.index[::-1], counts.values[::-1])
        ax.set_title(title)
        ax.set_xlabel("Features")
    fig.suptitle("Feature Mapping Flow (aggregated)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def write_markdown_summary(
    comparison: pd.DataFrame,
    *,
    out_path: Path,
    methodology: str,
    id_type_stats: dict | None = None,
    strategies: Sequence[str] | None = None,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Feature Mapping Summary (Legacy vs Enhancements)",
        "",
        "## Methodology",
        "",
        methodology,
        "",
        "## Definitions",
        "",
        "- **n_graph_nodes**: nodes of the mapped target type in the universal graph",
        "- **n_data_features**: feature columns in the omics table",
        "- **mapped_features / rate_from_data**: features with ≥1 in-graph hit / data features",
        "- **intersection_nodes / rate_to_network**: unique graph entities hit / graph nodes",
        "- **Enhancements** (BioMart / HGNC / UniProt / Overall): legacy first, then recover still-unmapped "
        "gene/transcript/protein features via that source profile. **Overall** unions BioMart+HGNC+UniProt.",
        "",
        "## Identifier detection",
        "",
        "```json",
        json.dumps(id_type_stats or {}, indent=2),
        "```",
        "",
        "## Coverage table (counts + both rates)",
        "",
    ]
    if comparison.empty:
        lines.append("_No comparison rows._")
    else:
        show = format_coverage_table(comparison, strategies=strategies)
        lines.append(show.to_markdown(index=False))
        delta_feat = comparison.get("delta_mapped_features", pd.Series([0])).fillna(0).sum()
        delta_nodes = comparison.get("delta_intersection_nodes", pd.Series([0])).fillna(0).sum()
        lines.extend(
            [
                "",
                "## Improvement highlights (baseline → overall/primary)",
                "",
                f"- Additional mapped features: **{int(delta_feat)}**",
                f"- Additional intersection graph nodes: **{int(delta_nodes)}**",
                "",
            ]
        )
        # Per-enhancement deltas vs legacy when available.
        per_src = []
        for s in ("biomart", "hgnc", "uniprot", "overall"):
            col = f"delta_{s}_mapped_features"
            if col in comparison.columns:
                per_src.append(
                    f"- `{s}` recovered features: **{int(comparison[col].fillna(0).sum())}**; "
                    f"intersection nodes: **{int(comparison.get(f'delta_{s}_intersection_nodes', pd.Series([0])).fillna(0).sum())}**"
                )
        if per_src:
            lines.extend(["## Recovery by enhancement source", "", *per_src, ""])
        lines.extend(
            [
                "## Recommendations",
                "",
                "- Prefer `mapping.strategy: optimized` with the **overall** profile for gene/transcript/protein.",
                "- Keep legacy for non-protein omics (miRNA, enhancer, promoter, metabolite).",
                "- ENST* inputs resolve via BioMart `ensembl_transcript_to_*.json` edge tables.",
                "",
            ]
        )
    out_path.write_text("\n".join(lines))


def write_html_summary(markdown_path: Path, html_path: Path, figure_paths: Iterable[Path]) -> None:
    """Lightweight HTML wrapper around the markdown summary + embedded figures."""
    html_path.parent.mkdir(parents=True, exist_ok=True)
    md = markdown_path.read_text() if markdown_path.exists() else ""
    # Minimal markdown→HTML (headings/tables already readable as pre).
    imgs = "\n".join(
        f'<img src="{p.name}" alt="{p.stem}" style="max-width:100%;margin:1em 0;"/>'
        for p in figure_paths
        if Path(p).exists()
    )
    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"/><title>Mapping Summary</title>
<style>
body {{ font-family: Georgia, serif; max-width: 1100px; margin: 2rem auto; padding: 0 1rem; }}
pre {{ white-space: pre-wrap; background: #f7f7f7; padding: 1rem; }}
table {{ border-collapse: collapse; }}
td, th {{ border: 1px solid #ccc; padding: 4px 8px; }}
</style></head><body>
<h1>Feature Mapping Summary</h1>
{imgs}
<pre>{md}</pre>
</body></html>
"""
    html_path.write_text(html)


def _omics_label(omics_file: str, node_type: str) -> str:
    mapping = {
        "transcript.csv": "Transcriptomics",
        "gene.csv": "Transcriptomics (gene)",
        "protein.csv": "Proteomics",
        "proteomics.csv": "Proteomics",
        "mirna.csv": "miRNA",
        "metabolite.csv": "Metabolomics",
    }
    return mapping.get(omics_file, node_type or omics_file)


def _base_uniprot(entity_id: str) -> str:
    eid = str(entity_id)
    for suffix in ("_transcript", "_protein"):
        if eid.endswith(suffix):
            return eid[: -len(suffix)]
    return eid
