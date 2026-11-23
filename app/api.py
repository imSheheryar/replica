"""Read-only FastAPI service for the Replica audit.

The API intentionally exposes only frozen artifacts.  It never fits a model or
accepts patient data, which keeps the public deployment free of PHI risk.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts"
TABLES = ARTIFACTS / "tables"
SITE = ROOT / "site"

app = FastAPI(
    title="Replica API",
    version="1.0.0",
    description="Which breast cancer subtypes survive external validation? A read-only audit API.",
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET"],
    allow_headers=["*"],
)


def _clean(value):
    if isinstance(value, dict):
        return {k: _clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_clean(v) for v in value]
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def load_results() -> dict:
    path = ARTIFACTS / "results.json"
    if not path.exists():
        raise HTTPException(status_code=503, detail="Analysis artifacts are not built yet")
    return _clean(json.loads(path.read_text()))


def load_table(name: str) -> list[dict]:
    path = TABLES / name
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"Missing table: {name}")
    return _clean(pd.read_csv(path).replace({float("nan"): None}).to_dict(orient="records"))


def compact_npc(npc: dict) -> dict:
    """Drop the large permutation-null arrays from API payloads."""
    out = {}
    for key, value in npc.items():
        if isinstance(value, dict):
            out[key] = {k: v for k, v in value.items() if k not in {"null", "null_combined"}}
        else:
            out[key] = value
    return out


@app.get("/health")
def health() -> dict:
    return {"status": "ok", "artifacts": (ARTIFACTS / "results.json").exists()}


@app.get("/api/summary")
def summary() -> dict:
    results = load_results()
    return {
        "meta": results["meta"],
        "decision": results["decision"],
        "cohorts": results["cohorts"],
        "npc": compact_npc(results["npc"]),
    }


@app.get("/api/scorecard")
def scorecard() -> dict:
    results = load_results()
    return {"scorecard": results["scorecard"], "primary_family": results["primary_family"]}


@app.get("/api/forest")
def forest() -> dict:
    return {"forest": load_table("forest.csv")}


@app.get("/api/km")
def km() -> dict:
    return {"curves": load_table("km_curves.csv")}


@app.get("/api/modules")
def modules() -> dict:
    results = load_results()
    return {"modules": results["modules"], "nmf": results["nmf"]}


@app.get("/api/npc")
def npc() -> dict:
    results = load_results()
    return {"npc": compact_npc(results["npc"]), "components": load_table("npc_components.csv")}


@app.get("/api/transport")
def transport() -> dict:
    results = load_results()
    return {"transport": results["transport"], "intclust": results["intclust"]}


@app.get("/api/methods")
def methods() -> dict:
    path = ROOT / "docs" / "methods.md"
    if not path.exists():
        raise HTTPException(status_code=404, detail="Methods document missing")
    return {"markdown": path.read_text()}


@app.get("/")
def index():
    path = SITE / "index.html"
    if not path.exists():
        raise HTTPException(status_code=503, detail="Dashboard not built yet")
    return FileResponse(path)


if (SITE / "assets").exists():
    app.mount("/assets", StaticFiles(directory=SITE / "assets"), name="assets")
