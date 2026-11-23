"""Nonparametric combination (NPC) and multiplicity control.

The central idea is to calibrate a combined test statistic against the *joint*
permutation distribution of its components.  That preserves the dependence
between endpoints and subtype contrasts without assuming a parametric copula.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from scipy.stats import norm, rankdata

Combiner = Literal["fisher", "liptak", "tippett"]


def _clip(p: np.ndarray, n: int) -> np.ndarray:
    return np.clip(p, 1.0 / (n + 1.0), 1.0 - 1.0 / (n + 1.0))


def combine_pvalues(pvalues: np.ndarray, method: Combiner = "liptak") -> float:
    """Combine marginal p-values into one statistic (larger = more evidence)."""
    p = np.asarray(pvalues, dtype=float)
    if p.size == 0:
        raise ValueError("At least one p-value is required")
    p = np.clip(p, 1e-12, 1 - 1e-12)
    if method == "fisher":
        return float(-2.0 * np.log(p).sum())
    if method == "liptak":
        return float(norm.ppf(1.0 - p).sum())
    if method == "tippett":
        return float(np.max(1.0 - p))
    raise ValueError(f"Unknown combiner: {method}")


def permutation_pvalue(observed: float, null: np.ndarray) -> float:
    """One-sided permutation p-value with the +1 correction."""
    null = np.asarray(null, dtype=float)
    return float((1.0 + np.sum(null >= observed)) / (null.size + 1.0))


def npc_test(
    observed: np.ndarray,
    null: np.ndarray,
    method: Combiner = "liptak",
) -> dict[str, object]:
    """Joint NPC test from a B x K matrix of permutation statistics.

    Parameters
    ----------
    observed:
        K observed component statistics.  Larger values mean stronger association.
    null:
        B x K statistics from the same permutation scheme.  Rows are permutations;
        columns are components.  The dependence among columns is preserved.
    method:
        Fisher, Liptak, or Tippett combining function.
    """
    observed = np.asarray(observed, dtype=float)
    null = np.asarray(null, dtype=float)
    if null.ndim != 2 or null.shape[1] != observed.size:
        raise ValueError("null must be B x K with K = len(observed)")
    n_perm = null.shape[0]

    component_p = np.array([permutation_pvalue(observed[j], null[:, j]) for j in range(observed.size)])
    # Rank-based marginal p-values for every permutation, preserving ties.
    component_p_null = np.empty_like(null)
    for j in range(observed.size):
        ranks = rankdata(null[:, j], method="average")
        component_p_null[:, j] = (n_perm + 1.0 - ranks) / (n_perm + 1.0)
    component_p_null = _clip(component_p_null, n_perm)

    observed_combined = combine_pvalues(component_p, method)
    null_combined = np.array([combine_pvalues(row, method) for row in component_p_null])
    combined_p = permutation_pvalue(observed_combined, null_combined)

    return {
        "method": method,
        "observed_components": observed.tolist(),
        "component_pvalues": component_p.tolist(),
        "observed_combined": observed_combined,
        "combined_pvalue": combined_p,
        "null_combined": null_combined,
        "n_permutations": int(n_perm),
    }


def bh(pvalues: np.ndarray, q: float = 0.05) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def by(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Yekutieli adjusted p-values (arbitrary dependence)."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    harmonic = np.sum(1.0 / np.arange(1, n + 1))
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n * harmonic / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def maxt_adjusted(observed: np.ndarray, null: np.ndarray) -> np.ndarray:
    """Step-down-free maxT adjusted p-values from a B x K permutation matrix."""
    observed = np.asarray(observed, dtype=float)
    null = np.asarray(null, dtype=float)
    maxima = null.max(axis=1)
    return np.array([(1.0 + np.sum(maxima >= value)) / (null.shape[0] + 1.0) for value in observed])


def fdr_mask(pvalues: np.ndarray, q: float = 0.05, method: str = "bh") -> np.ndarray:
    adjusted = bh(pvalues, q) if method == "bh" else by(pvalues)
    return adjusted <= q
