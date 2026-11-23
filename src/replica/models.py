"""Expression modules: NMF in the discovery cohort, locked scoring externally."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score


@dataclass
class NMFResult:
    rank: int
    model: NMF
    W: pd.DataFrame
    H: pd.DataFrame
    reconstruction_error: float
    silhouette: float


def select_variable_genes(expr: pd.DataFrame, n: int = 500) -> list[str]:
    variance = expr.var(axis=0).sort_values(ascending=False)
    return variance.head(n).index.tolist()


def _nonnegative(matrix: np.ndarray) -> np.ndarray:
    shifted = matrix - np.nanmin(matrix)
    return np.clip(shifted, 0, None) + 1e-9


def fit_nmf(
    expr: pd.DataFrame,
    ranks: list[int] | None = None,
    max_iter: int = 600,
    seed: int = 20261123,
) -> tuple[NMFResult, pd.DataFrame]:
    """Fit NMF across candidate ranks and select by silhouette stability."""
    ranks = ranks or [3, 4, 5, 6, 7, 8]
    genes = select_variable_genes(expr, min(500, expr.shape[1]))
    x = _nonnegative(expr[genes].to_numpy(dtype=float))
    results: list[NMFResult] = []
    for rank in ranks:
        model = NMF(
            n_components=rank,
            init="nndsvda",
            random_state=seed,
            max_iter=max_iter,
            beta_loss="frobenius",
            solver="mu",
        )
        w = model.fit_transform(x)
        h = model.components_
        labels = w.argmax(axis=1)
        sil = silhouette_score(w, labels) if len(set(labels)) > 1 else -1.0
        results.append(
            NMFResult(
                rank=rank,
                model=model,
                W=pd.DataFrame(w, index=expr.index, columns=[f"M{i+1}" for i in range(rank)]),
                H=pd.DataFrame(h, columns=genes, index=[f"M{i+1}" for i in range(rank)]),
                reconstruction_error=float(model.reconstruction_err_),
                silhouette=float(sil),
            )
        )
    best_sil = max(r.silhouette for r in results)
    # Prefer the smallest rank within 0.02 silhouette of the best.
    eligible = [r for r in results if r.silhouette >= best_sil - 0.02]
    chosen = min(eligible, key=lambda r: r.rank)
    comparison = pd.DataFrame(
        [
            {
                "rank": r.rank,
                "reconstruction_error": r.reconstruction_error,
                "silhouette": r.silhouette,
                "selected": r.rank == chosen.rank,
            }
            for r in results
        ]
    )
    return chosen, comparison


def top_module_genes(result: NMFResult, n: int = 30) -> dict[str, list[str]]:
    modules: dict[str, list[str]] = {}
    for module in result.H.index:
        modules[module] = result.H.loc[module].sort_values(ascending=False).head(n).index.tolist()
    return modules


def module_scores(expr: pd.DataFrame, modules: dict[str, list[str]]) -> pd.DataFrame:
    scores = {}
    for module, genes in modules.items():
        available = [g for g in genes if g in expr.columns]
        if available:
            scores[module] = expr[available].mean(axis=1)
    return pd.DataFrame(scores)
