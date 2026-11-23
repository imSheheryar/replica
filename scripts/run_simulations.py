"""Simulation study: type I error, power, and FDR under dependent censored endpoints.

The simulation is deliberately small enough to run in CI-adjacent environments.
It is a methods check, not a breast cancer result.
"""

from __future__ import annotations


import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

from replica.config import PATHS
from replica.stats import bh, npc_test
from replica.survival import logrank_statistic

SEED = 20261123


def correlated_times(rng, n, p, rho, shift, treatment, censor_scale=60.0):
    """Multivariate t errors -> correlated exponential survival times."""
    cov = np.full((p, p), rho)
    np.fill_diagonal(cov, 1.0)
    chol = np.linalg.cholesky(cov)
    z = rng.standard_t(df=5, size=(n, p)) @ chol.T
    u = np.clip(norm.cdf(z), 1e-6, 1 - 1e-6)
    hazard = 0.06 * np.exp(shift * treatment[:, None])
    event_time = -np.log(u) / hazard
    censor = rng.exponential(censor_scale, size=(n, p))
    time = np.minimum(event_time, censor)
    event = (event_time <= censor).astype(int)
    return time, event


def replicate(rng, n=160, p=5, rho=0.45, shift=0.0, n_perm=100):
    treatment = rng.integers(0, 2, n)
    time, event = correlated_times(rng, n, p, rho, shift, treatment)
    observed = np.array([logrank_statistic(time[:, j], event[:, j], treatment) for j in range(p)])
    null = np.empty((n_perm, p))
    for b in range(n_perm):
        perm = rng.permutation(treatment)
        for j in range(p):
            null[b, j] = logrank_statistic(time[:, j], event[:, j], perm)
    npc = npc_test(observed, null, method="fisher")
    marginal_p = chi2.sf(observed, 1)
    return {
        "npc_p": npc["combined_pvalue"],
        "marginal_bh_p": float(bh(marginal_p).min()),
        "bonferroni_p": float(min(1.0, marginal_p.min() * p)),
        "min_marginal_p": float(marginal_p.min()),
    }


def run_scenario(shift: float, n_rep: int = 150, n_perm: int = 80) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + int(shift * 1000))
    rows = []
    for rep in range(n_rep):
        result = replicate(rng, shift=shift, n_perm=n_perm)
        result["replicate"] = rep
        result["shift"] = shift
        rows.append(result)
    return pd.DataFrame(rows)


def fdr_scenario(n_features: int = 80, n_signal: int = 16, n_rep: int = 40, shift: float = 0.45) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + 99)
    rows = []
    for rep in range(n_rep):
        treatment = rng.integers(0, 2, 160)
        pvals_npc = []
        pvals_marginal = []
        is_signal = np.zeros(n_features, dtype=bool)
        is_signal[:n_signal] = True
        rng.shuffle(is_signal)
        for f in range(n_features):
            time, event = correlated_times(
                rng, 160, 5, 0.45, shift if is_signal[f] else 0.0, treatment
            )
            stats = np.array([logrank_statistic(time[:, j], event[:, j], treatment) for j in range(5)])
            p = chi2.sf(stats, 1)
            pvals_marginal.append(float(bh(p).min()))
            # Asymptotic Fisher combination; the permutation calibration is
            # validated separately in tests and the main analysis.
            fisher = float(chi2.sf(-2 * np.log(np.clip(p, 1e-12, 1)).sum(), 10))
            pvals_npc.append(fisher)
        for name, pvals in [("marginal_BH", pvals_marginal), ("NPC_Fisher", pvals_npc)]:
            q = bh(np.asarray(pvals), q=0.05)
            discoveries = q <= 0.05
            false = discoveries & ~is_signal
            rows.append(
                {
                    "replicate": rep,
                    "method": name,
                    "discoveries": int(discoveries.sum()),
                    "false_discoveries": int(false.sum()),
                    "fdp": float(false.sum() / max(1, discoveries.sum())),
                    "power": float((discoveries & is_signal).sum() / max(1, is_signal.sum())),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    PATHS.ensure()
    rows = []
    for shift in [0.0, 0.35, 0.55]:
        frame = run_scenario(shift)
        rows.append(frame)
        print(f"shift={shift}: npc rejection={(frame.npc_p < 0.05).mean():.3f}")
    sim = pd.concat(rows, ignore_index=True)
    sim.to_csv(PATHS.tables / "simulation_replicates.csv", index=False)

    summary = (
        sim.groupby("shift")
        .agg(
            n=("replicate", "count"),
            npc_type1_or_power=("npc_p", lambda x: float((x < 0.05).mean())),
            marginal_bh_rate=("marginal_bh_p", lambda x: float((x < 0.05).mean())),
            bonferroni_rate=("bonferroni_p", lambda x: float((x < 0.05).mean())),
        )
        .reset_index()
    )
    summary.to_csv(PATHS.tables / "simulation_summary.csv", index=False)

    fdr = fdr_scenario()
    fdr.to_csv(PATHS.tables / "simulation_fdr_replicates.csv", index=False)
    fdr_summary = (
        fdr.groupby("method")
        .agg(
            mean_discoveries=("discoveries", "mean"),
            mean_fdp=("fdp", "mean"),
            mean_power=("power", "mean"),
        )
        .reset_index()
    )
    fdr_summary.to_csv(PATHS.tables / "simulation_fdr_summary.csv", index=False)
    print(summary.to_string(index=False))
    print(fdr_summary.to_string(index=False))


if __name__ == "__main__":
    main()
