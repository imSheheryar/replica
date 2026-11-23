from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "config" / "analysis.yaml"


@dataclass(frozen=True)
class Paths:
    root: Path = ROOT
    raw: Path = ROOT / "data" / "raw"
    interim: Path = ROOT / "data" / "interim"
    processed: Path = ROOT / "data" / "processed"
    artifacts: Path = ROOT / "artifacts"
    figures: Path = ROOT / "artifacts" / "figures"
    tables: Path = ROOT / "artifacts" / "tables"
    models: Path = ROOT / "artifacts" / "models"
    site: Path = ROOT / "site"

    def ensure(self) -> None:
        for path in (
            self.raw,
            self.interim,
            self.processed,
            self.artifacts,
            self.figures,
            self.tables,
            self.models,
            self.site,
        ):
            path.mkdir(parents=True, exist_ok=True)


def load_config(path: Path | None = None) -> dict[str, Any]:
    with (path or CONFIG_PATH).open() as handle:
        return yaml.safe_load(handle)


PATHS = Paths()
