"""Path configuration for universal graph and offline ID-mapping tables."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


@dataclass
class MapperPaths:
    """Filesystem locations required by the mapper."""

    universal_graph: Path
    biomart_dir: Path
    hgnc_dir: Path
    uniprot_dir: Path

    def validate(self) -> None:
        missing = []
        if not self.universal_graph.is_file():
            missing.append(f"universal_graph: {self.universal_graph}")
        for name, p in (
            ("biomart_dir", self.biomart_dir),
            ("hgnc_dir", self.hgnc_dir),
            ("uniprot_dir", self.uniprot_dir),
        ):
            if not p.is_dir():
                missing.append(f"{name}: {p}")
        if missing:
            raise FileNotFoundError(
                "Missing mapper inputs. Edit config/paths.yaml or set UGM_* env vars.\n"
                + "\n".join(f"  - {m}" for m in missing)
            )


def _env_path(key: str) -> Path | None:
    val = os.environ.get(key, "").strip()
    return Path(val).expanduser() if val else None


def load_paths(config_path: str | Path | None = None) -> MapperPaths:
    """Load paths from YAML (default: package config/paths.yaml) with UGM_* env overrides."""
    pkg_root = Path(__file__).resolve().parent.parent
    cfg_file = Path(config_path) if config_path else pkg_root / "config" / "paths.yaml"
    data: dict[str, Any] = {}
    if cfg_file.is_file() and yaml is not None:
        data = yaml.safe_load(cfg_file.read_text()) or {}

    def pick(key: str, env: str, default: str = "") -> Path:
        override = _env_path(env)
        if override is not None:
            return override
        raw = data.get(key, default)
        p = Path(str(raw)).expanduser()
        return p if p.is_absolute() else (pkg_root / p).resolve()

    return MapperPaths(
        universal_graph=pick("universal_graph", "UGM_UNIVERSAL_GRAPH"),
        biomart_dir=pick("biomart_dir", "UGM_BIOMART_DIR"),
        hgnc_dir=pick("hgnc_dir", "UGM_HGNC_DIR"),
        uniprot_dir=pick("uniprot_dir", "UGM_UNIPROT_DIR"),
    )
