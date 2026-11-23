"""Build the static Replica evidence dashboard from frozen artifacts."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts"
TABLES = ARTIFACTS / "tables"
SITE = ROOT / "site"
ASSETS = SITE / "assets"

INK = "#12212b"
MUTED = "#5f6f78"
PAPER = "#fbfaf7"
TEAL = "#1f7a72"
CORAL = "#d45d4c"
GOLD = "#c99a3c"
BLUE = "#3b6ea8"
GRID = "#e6e2d9"

pio.templates.default = "plotly_white"
FONT = "Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif"


def style(fig, height=420, margin=None):
    fig.update_layout(
        height=height,
        margin=margin or dict(l=50, r=30, t=50, b=50),
        font=dict(family=FONT, color=INK, size=12),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        hoverlabel=dict(font_family=FONT),
    )
    fig.update_xaxes(gridcolor=GRID, zeroline=False)
    fig.update_yaxes(gridcolor=GRID, zeroline=False)
    return fig


def fig_html(fig, div_id):
    return pio.to_html(
        fig,
        full_html=False,
        include_plotlyjs=False,
        div_id=div_id,
        config={"displayModeBar": False, "responsive": True},
    )


def forest_figure(forest: pd.DataFrame):
    data = forest[forest["term"].isin(["Basal", "Her2", "LumB"])].copy()
    data["label"] = data["cohort"].str.upper() + " · " + data["endpoint"]
    colors = {"tcga": TEAL, "metabric": CORAL}
    fig = go.Figure()
    for term, symbol in [("Basal", "circle"), ("Her2", "diamond"), ("LumB", "square")]:
        subset = data[data["term"] == term]
        for cohort in ["tcga", "metabric"]:
            s = subset[subset["cohort"] == cohort]
            if s.empty:
                continue
            fig.add_trace(
                go.Scatter(
                    x=s["hr"],
                    y=s["label"],
                    mode="markers",
                    name=f"{term} ({cohort.upper()})",
                    marker=dict(size=9, symbol=symbol, color=colors[cohort], line=dict(width=1, color="white")),
                    error_x=dict(
                        type="data",
                        symmetric=False,
                        array=s["upper"] - s["hr"],
                        arrayminus=s["hr"] - s["lower"],
                        color=colors[cohort],
                        thickness=1.2,
                    ),
                    hovertemplate="%{y}<br>HR %{x:.2f}<extra></extra>",
                )
            )
    fig.add_vline(x=1, line_dash="dash", line_color=MUTED, line_width=1)
    fig.update_xaxes(type="log", title="Hazard ratio vs LumA (log scale)")
    fig.update_yaxes(title="")
    fig.update_layout(legend=dict(orientation="h", y=-0.18, x=0))
    return style(fig, height=430)


def km_figure(km: pd.DataFrame):
    specs = [("tcga", "OS"), ("tcga", "DFS"), ("metabric", "OS"), ("metabric", "RFS")]
    titles = ["TCGA · overall survival", "TCGA · disease-free survival", "METABRIC · overall survival", "METABRIC · relapse-free survival"]
    fig = make_subplots(rows=2, cols=2, subplot_titles=titles, shared_xaxes=False, vertical_spacing=0.16, horizontal_spacing=0.1)
    palette = {"LumA": TEAL, "LumB": GOLD, "Basal": CORAL, "Her2": BLUE}
    for i, (cohort, endpoint) in enumerate(specs):
        row, col = i // 2 + 1, i % 2 + 1
        subset = km[(km["cohort"] == cohort) & (km["endpoint"] == endpoint)]
        for group in ["LumA", "LumB", "Basal", "Her2"]:
            g = subset[subset["group"] == group]
            if g.empty:
                continue
            fig.add_trace(
                go.Scatter(
                    x=g["time"],
                    y=g["survival"],
                    mode="lines",
                    name=group,
                    line=dict(color=palette[group], width=2),
                    legendgroup=group,
                    showlegend=(i == 0),
                    hovertemplate=f"{group}<br>%{{x:.0f}} months<br>%{{y:.0%}}<extra></extra>",
                ),
                row=row,
                col=col,
            )
    fig.update_xaxes(title_text="Months", row=2, col=1)
    fig.update_xaxes(title_text="Months", row=2, col=2)
    fig.update_yaxes(title_text="Survival", range=[0, 1.02], row=1, col=1)
    fig.update_yaxes(title_text="Survival", range=[0, 1.02], row=2, col=1)
    fig.update_layout(legend=dict(orientation="h", y=-0.12, x=0))
    return style(fig, height=650)


def incremental_figure(inc: pd.DataFrame):
    data = inc[inc["endpoint"] == "OS"].copy()
    labels = ["TCGA", "METABRIC"]
    base = data[data["cohort"] == "tcga"]["cindex_base"].iloc[0]
    full = data[data["cohort"] == "tcga"]["cindex_full"].iloc[0]
    base_m = data[data["cohort"] == "metabric"]["cindex_base"].iloc[0]
    full_m = data[data["cohort"] == "metabric"]["cindex_full"].iloc[0]
    fig = go.Figure()
    fig.add_trace(go.Bar(name="Clinical", x=labels, y=[base, base_m], marker_color="#b9c7c9",
                         error_y=dict(type="data", array=[data.iloc[0]["cindex_base_upper"] - base, data.iloc[1]["cindex_base_upper"] - base_m], color=MUTED)))
    fig.add_trace(go.Bar(name="Clinical + PAM50", x=labels, y=[full, full_m], marker_color=TEAL,
                         error_y=dict(type="data", array=[data.iloc[0]["cindex_full_upper"] - full, data.iloc[1]["cindex_full_upper"] - full_m], color=TEAL)))
    fig.add_hline(y=0.5, line_dash="dot", line_color=MUTED, line_width=1)
    fig.add_annotation(x=1.02, y=0.5, xref="paper", text="chance", showarrow=False, font=dict(size=10, color=MUTED), xanchor="left")
    fig.update_yaxes(title="Out-of-fold C-index (axis starts at chance)", range=[0.5, 0.85])
    fig.update_layout(barmode="group", legend=dict(orientation="h", y=1.08, x=0))
    for i, (b, f) in enumerate([(base, full), (base_m, full_m)]):
        fig.add_annotation(x=labels[i], y=max(f, b) + 0.035, text=f"ΔC {f-b:+.3f}", showarrow=False, font=dict(size=11, color=INK))
    return style(fig, height=390)


def calibration_figure(cal: pd.DataFrame):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(color=MUTED, dash="dash"), name="Ideal"))
    fig.add_trace(
        go.Scatter(
            x=cal["predicted_survival"],
            y=cal["observed_survival"],
            mode="markers+lines",
            marker=dict(size=11, color=CORAL, line=dict(width=1, color="white")),
            line=dict(color=CORAL, width=2),
            name="Locked TCGA model",
            text=[f"Q{int(q)}" for q in cal["quartile"]],
            textposition="top center",
            textfont=dict(size=10, color=CORAL),
            hovertemplate="Predicted %{x:.0%}<br>Observed %{y:.0%}<extra></extra>",
        )
    )
    fig.update_xaxes(title="Predicted 60-month survival (TCGA model)", tickformat=".0%", range=[0, 1])
    fig.update_yaxes(title="Observed 60-month survival (METABRIC)", tickformat=".0%", range=[0, 1])
    fig.update_layout(legend=dict(orientation="h", y=1.08, x=0))
    return style(fig, height=390)


def npc_figure(npc: dict):
    null = np.asarray(npc["null_combined"], dtype=float)
    observed = npc["observed_combined"]
    p = npc["combined_pvalue"]
    fig = go.Figure()
    fig.add_trace(go.Histogram(x=null, nbinsx=50, marker_color="rgba(31,122,114,0.65)", name="Permutation null"))
    fig.add_vline(x=observed, line_color=CORAL, line_width=2, line_dash="dash")
    p_text = f"p ≤ {p:.4f}" if p <= 1 / (npc.get("n_permutations", 2000) + 1) else f"p = {p:.4f}"
    fig.add_annotation(x=observed, y=1, yref="paper", text=f"observed<br>{p_text}", showarrow=False, xanchor="left", font=dict(color=CORAL, size=11))
    fig.update_xaxes(title="Liptak combined statistic (4 components)")
    fig.update_yaxes(title="Permutations")
    return style(fig, height=370)


def module_heatmap(modules: pd.DataFrame):
    data = modules.copy()
    data["module"] = data["module"].astype(str)
    cohorts = ["tcga", "metabric"]
    fig = make_subplots(rows=1, cols=2, subplot_titles=["TCGA", "METABRIC"], horizontal_spacing=0.08)
    for i, cohort in enumerate(cohorts):
        subset = data[data["cohort"] == cohort]
        if subset.empty:
            continue
        # The module table stores one row per module; subtype association is a
        # single p-value, so the heatmap shows signed effect direction and -log10 p.
        z = subset["cox_hr_per_sd"].to_numpy().reshape(-1, 1)
        fig.add_trace(
            go.Heatmap(
                z=z,
                y=subset["module"],
                x=["HR per SD"],
                colorscale=[[0, "#7fa7c9"], [0.5, "#f4f1ea"], [1, CORAL]],
                zmid=1,
                zmin=0.8,
                zmax=1.2,
                text=[[f"{v:.2f}"] for v in z[:, 0]],
                texttemplate="%{text}",
                textfont=dict(size=14, color=INK),
                colorbar=dict(title="HR", len=0.7, x=0.46 if i == 0 else 1.02),
                hovertemplate="%{y}<br>HR per SD %{z:.2f}<extra></extra>",
                showscale=(i == 1),
            ),
            row=1,
            col=i + 1,
        )
    fig.update_xaxes(showticklabels=False)
    return style(fig, height=360, margin=dict(l=60, r=60, t=50, b=30))


def rank_figure(nmf: dict):
    comp = pd.DataFrame(nmf["comparison"])
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=comp["rank"], y=comp["silhouette"], mode="lines+markers", line=dict(color=TEAL, width=2), marker=dict(size=9), name="Silhouette"))
    selected = comp[comp["selected"].astype(bool)]
    if not selected.empty:
        fig.add_trace(go.Scatter(x=selected["rank"], y=selected["silhouette"], mode="markers", marker=dict(size=15, color=CORAL, symbol="circle-open", line=dict(width=2)), name="Selected"))
    fig.update_xaxes(title="NMF rank (k)", dtick=1)
    fig.update_yaxes(title="Mean silhouette of W")
    fig.update_layout(legend=dict(orientation="h", y=1.08, x=0))
    return style(fig, height=330)


def scorecard_table(scorecard: pd.DataFrame) -> str:
    rows = []
    for _, r in scorecard.iterrows():
        p = f"{r['p_value']:.2e}" if pd.notna(r["p_value"]) else "—"
        q = f"{r['p_bh']:.3f}" if pd.notna(r.get("p_bh", np.nan)) else "—"
        rows.append(
            f"<tr><td><strong>{r['cohort'].upper()}</strong></td><td>{r['endpoint']}</td>"
            f"<td>{int(r['n']):,}</td><td>{int(r['events']):,}</td><td>{r['cindex']:.3f}</td>"
            f"<td>{p}</td><td>{q}</td></tr>"
        )
    return "".join(rows)


def module_table(modules: pd.DataFrame) -> str:
    rows = []
    for _, r in modules.iterrows():
        q = f"{r['cox_q']:.3f}" if pd.notna(r["cox_q"]) else "—"
        rows.append(
            f"<tr><td><strong>{r['module']}</strong></td><td>{r['cohort'].upper()}</td>"
            f"<td>{r['cox_hr_per_sd']:.2f}</td><td>{r['cox_p']:.2e}</td><td>{q}</td></tr>"
        )
    return "".join(rows)


def sensitivity_table(rows: list[dict]) -> str:
    if not rows:
        return ""
    body = []
    for row in rows:
        body.append(
            f"<tr><td><strong>{row['cohort'].upper()}</strong></td><td>{row['endpoint']}</td>"
            f"<td>{int(row['n']):,}</td><td>{int(row['events']):,}</td>"
            f"<td>{row['cindex']:.3f}</td><td>{row['p_value']:.2e}</td></tr>"
        )
    return "".join(body)


def simulation_figure(summary: pd.DataFrame):
    fig = go.Figure()
    series = [
        ("npc_type1_or_power", "NPC (Fisher)", TEAL),
        ("marginal_bh_rate", "Marginal BH", "#8a9a5b"),
        ("bonferroni_rate", "Bonferroni", CORAL),
    ]
    for column, label, color in series:
        fig.add_trace(
            go.Scatter(
                x=summary["shift"],
                y=summary[column],
                mode="lines+markers",
                name=label,
                line=dict(color=color, width=2),
                marker=dict(size=8),
                hovertemplate="shift %{x:.2f}<br>rejection %{y:.1%}<extra></extra>",
            )
        )
    fig.add_hline(y=0.05, line_dash="dot", line_color=MUTED, line_width=1)
    fig.add_annotation(x=0.02, y=0.05, text="nominal 0.05", showarrow=False, font=dict(size=10, color=MUTED), xanchor="left", yshift=8)
    fig.update_xaxes(title="Log-hazard shift (0 = null)", tickmode="array", tickvals=summary["shift"].tolist(), ticktext=["null" if s == 0 else f"{s:.2f}" for s in summary["shift"]])
    fig.update_yaxes(title="Rejection rate at α = 0.05", range=[0, 1], tickformat=".0%")
    fig.update_layout(legend=dict(orientation="h", y=1.12, x=0))
    return style(fig, height=380)


def simulation_table(fdr: pd.DataFrame) -> str:
    rows = []
    for _, r in fdr.iterrows():
        method = "NPC (Fisher)" if r["method"] == "NPC_Fisher" else r["method"].replace("_", " ")
        rows.append(
            f"<tr><td><strong>{method}</strong></td><td>{r['mean_discoveries']:.1f}</td>"
            f"<td>{r['mean_fdp']:.3f}</td><td>{r['mean_power']:.1%}</td></tr>"
        )
    return "".join(rows)


def build() -> Path:
    ASSETS.mkdir(parents=True, exist_ok=True)
    results = json.loads((ARTIFACTS / "results.json").read_text())
    scorecard = pd.read_csv(TABLES / "scorecard.csv")
    forest = pd.read_csv(TABLES / "forest.csv")
    km = pd.read_csv(TABLES / "km_curves.csv")
    inc = pd.read_csv(TABLES / "incremental.csv")
    cal = pd.read_csv(TABLES / "calibration.csv")
    modules = pd.read_csv(TABLES / "modules.csv")
    sim_path = TABLES / "simulation_summary.csv"
    fdr_path = TABLES / "simulation_fdr_summary.csv"
    sim_summary = pd.read_csv(sim_path) if sim_path.exists() else None
    fdr_summary = pd.read_csv(fdr_path) if fdr_path.exists() else None
    sensitivity = results.get("sensitivity", [])
    nmf = results["nmf"]
    decision = results["decision"]
    cohorts = results["cohorts"]

    figs = {
        "forest": fig_html(forest_figure(forest), "fig-forest"),
        "km": fig_html(km_figure(km), "fig-km"),
        "incremental": fig_html(incremental_figure(inc), "fig-incremental"),
        "calibration": fig_html(calibration_figure(cal), "fig-calibration"),
        "npc": fig_html(npc_figure(results["npc"]["global"]), "fig-npc"),
        "modules": fig_html(module_heatmap(modules), "fig-modules"),
        "rank": fig_html(rank_figure(nmf), "fig-rank"),
    }
    if sim_summary is not None:
        figs["simulation"] = fig_html(simulation_figure(sim_summary), "fig-simulation")

    total_n = cohorts["tcga"]["n_raw"] + cohorts["metabric"]["n_raw"]
    os_n = cohorts["tcga"]["n_analysis_os"] + cohorts["metabric"]["n_analysis_os"]
    badge_class = "good" if decision["replicated"] else "mixed"
    badge_text = decision["label"]
    global_p = results["npc"]["global"]["combined_pvalue"]
    global_components = results["npc"]["global"]["component_pvalues"]
    transport = results["transport"]
    delta = transport["table"][1]["delta_cindex"] if len(transport["table"]) > 1 else float("nan")
    delta_lo = transport["table"][1]["delta_lower"] if len(transport["table"]) > 1 else float("nan")
    delta_hi = transport["table"][1]["delta_upper"] if len(transport["table"]) > 1 else float("nan")
    c_base = transport["table"][0]["cindex"] if transport["table"] else float("nan")
    c_full = transport["table"][1]["cindex"] if len(transport["table"]) > 1 else float("nan")
    intclust = results["intclust"]
    intclust_text = (
        f"METABRIC internal benchmark: clinical {intclust[0]['cindex']:.3f} → +PAM50 "
        f"{intclust[1]['cindex']:.3f} → +PAM50+IntClust {intclust[2]['cindex']:.3f}; "
        f"IntClust adds {intclust[2]['delta_vs_pam50']:+.4f} over PAM50. This is internal, not external validation."
    )
    global_p_text = f"p ≤ {global_p:.4f}" if global_p <= 1 / (results["npc"]["global"].get("n_permutations", 2000) + 1) else f"p = {global_p:.4f}"
    component_text = ", ".join(f"{p:.3f}" if p > 0.001 else "≤0.001" for p in global_components)
    if decision["replicated"]:
        decision_detail = "The replication rule was met across both OS cohorts."
    elif decision["label"] == "Cohort-specific":
        decision_detail = "The replication rule was not met across both OS cohorts: one cohort separated and the other did not."
    else:
        decision_detail = "The replication rule was not met across both OS cohorts."

    narrative = """
    <p class="lede">Breast cancer intrinsic subtypes are often presented as portable biology. Replica tests that claim the way an external audit should: freeze the published calls, lock the model, move it to an independent cohort, and ask whether it still adds information beyond age, stage, and receptor status.</p>
    <p>The audit compares the published PAM50/Claudin-low calls in <strong>TCGA-BRCA</strong> (reference cohort) and <strong>METABRIC</strong> (external validation). It reports the subtype association, the incremental value over a clinical baseline, and a censoring-aware nonparametric combination across OS and DFS/RFS. IntClust is used as an internal METABRIC benchmark because no external IntClust cohort exists.</p>
    """

    findings = f"""
    <div class="finding"><span class="finding-kicker">Decision</span><h3>{badge_text}</h3><p>{decision_detail} Incremental value {'was present' if decision['incremental'] else 'was not established'} under the prespecified ΔC rule.</p></div>
    <div class="finding"><span class="finding-kicker">Global NPC</span><h3>{global_p_text}</h3><p>Liptak combination of four dependent components: TCGA OS, TCGA DFS, METABRIC OS, METABRIC RFS. Components: {component_text}. The two METABRIC components are at the floor; the TCGA components are not. This is a global association test, not a replication test.</p></div>
    <div class="finding"><span class="finding-kicker">Transport</span><h3>ΔC {delta:+.3f}</h3><p>Locked TCGA clinical + PAM50 model applied to METABRIC. C-index {c_base:.3f} → {c_full:.3f}; paired bootstrap 95% CI [{delta_lo:+.4f}, {delta_hi:+.4f}]. Calibration slope {transport['calibration_slope']:.2f}.</p></div>
    <div class="finding"><span class="finding-kicker">Modules</span><h3>k = {nmf['selected_rank']}</h3><p>NMF modules derived in TCGA from a pre-specified 995-gene Hallmark panel and scored in METABRIC without using METABRIC outcomes.</p></div>
    """
    caveat = decision.get("caveat", "")
    if sim_summary is not None and fdr_summary is not None:
        simulation_section = f"""
  <section>
    <div class="section-head"><span class="section-no">10</span><h2>Methods check: operating characteristics</h2><p class="section-note">A simulation, not a breast cancer result. At shift 0 the rejection rate estimates Type I error; positive shifts estimate power. The simulation uses Fisher combination for speed; the cohort-level primary is permutation-calibrated Liptak.</p></div>
    <div class="grid-2">
      <div class="figure">{figs['simulation']}</div>
      <div class="panel"><h3>FDR under dependent features</h3>
        <div class="table-wrap"><table><thead><tr><th>Method</th><th>Discoveries</th><th>Mean FDP</th><th>Power</th></tr></thead><tbody>{simulation_table(fdr_summary)}</tbody></table></div>
        <p class="small">The asymptotic Fisher variant is anti-conservative under dependent endpoints (mean FDP 0.108), which is exactly why the cohort analysis calibrates Liptak against joint permutations instead of using asymptotic combination. The marginal BH arm tests one endpoint at a time and is not directly comparable to the multi-endpoint NPC. This is a methods check, not a second cohort p-value.</p>
      </div>
    </div>
  </section>
"""
        limits_no, repro_no = "11", "12"
    else:
        simulation_section = ""
        limits_no, repro_no = "10", "11"

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Replica — Which breast cancer subtypes survive external validation?</title>
<meta name="description" content="A replication audit of breast cancer subtype systems across TCGA-BRCA and METABRIC, with censoring-aware NPC, FDR control, and external validation.">
<link rel="icon" type="image/svg+xml" href="assets/favicon.svg">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=DM+Serif+Display&family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<link rel="stylesheet" href="assets/style.css">
<script src="assets/plotly.min.js"></script>
</head>
<body>
<header class="hero">
  <div class="wrap">
    <div class="topline"><span class="mark">R / 01</span><span>Replication audit · biomedical data science</span><span class="topline-right">TCGA-BRCA × METABRIC</span></div>
    <p class="eyebrow">Which breast cancer subtypes survive external validation?</p>
    <h1>Replica<span class="dot">.</span></h1>
    <p class="subtitle">A pre-specified, internally frozen audit of PAM50, clinical receptors, and expression modules across two independent breast cancer cohorts — with a censoring-aware nonparametric combination framework.</p>
    <div class="hero-meta">
      <span class="pill {badge_class}">{badge_text}</span>
      <span class="pill ghost">{total_n:,} patients</span>
      <span class="pill ghost">{os_n:,} in OS models</span>
      <span class="pill ghost">{results['meta']['n_permutations_primary']:,} permutations</span>
    </div>
    <div class="hero-grid">
      <div><span class="stat">{cohorts['tcga']['n_raw']:,}</span><span class="stat-label">TCGA-BRCA samples</span></div>
      <div><span class="stat">{cohorts['metabric']['n_raw']:,}</span><span class="stat-label">METABRIC patients</span></div>
      <div><span class="stat">{len(scorecard)}</span><span class="stat-label">cohort × endpoint tests</span></div>
      <div><span class="stat">{nmf['selected_rank']}</span><span class="stat-label">locked expression modules</span></div>
    </div>
  </div>
</header>

<main class="wrap">
  <section class="intro">
    <div class="section-head"><span class="section-no">01</span><h2>The question</h2></div>
    {narrative}
  </section>

  <section>
    <div class="section-head"><span class="section-no">02</span><h2>Headline evidence</h2></div>
    <div class="findings">{findings}</div>
  </section>

  <section>
    <div class="section-head"><span class="section-no">03</span><h2>Replication scorecard</h2><p class="section-note">Primary family uses Benjamini–Hochberg at q = 0.05 across the four prespecified OS tests. C-index is unadjusted within cohort; the external transport estimate is shown separately. A p of 0.0005 is the Monte Carlo floor at 2,000 permutations (1/2001), not an exact population p-value.</p></div>
    <div class="table-wrap"><table><thead><tr><th>Cohort</th><th>Endpoint</th><th>n</th><th>Events</th><th>C-index</th><th>Log-rank p</th><th>BH q</th></tr></thead><tbody>{scorecard_table(scorecard)}</tbody></table></div>
    <div class="callout"><strong>What the null result does and does not mean.</strong> {caveat}</div>
  </section>

  <section>
    <div class="section-head"><span class="section-no">03b</span><h2>Sensitivity: the complete-case filter matters</h2><p class="section-note">The primary cohort requires age, stage, ER, and PR. This sensitivity uses every patient with a four-class subtype and a valid endpoint, without the clinical complete-case filter. It is not a rescue; it shows how much the exclusion changes the answer.</p></div>
    <div class="table-wrap"><table><thead><tr><th>Cohort</th><th>Endpoint</th><th>n</th><th>Events</th><th>C-index</th><th>Permutation p</th></tr></thead><tbody>{sensitivity_table(sensitivity)}</tbody></table></div>
  </section>

  <section>
    <div class="section-head"><span class="section-no">04</span><h2>Does subtype add anything?</h2><p class="section-note">Out-of-fold C-index from 5-fold cross-validation within each cohort. ΔC is paired; the bootstrap interval is the primary incremental estimand.</p></div>
    <div class="figure wide">{figs['incremental']}</div>
  </section>

  <section>
    <div class="section-head"><span class="section-no">05</span><h2>Hazard ratios by subtype</h2><p class="section-note">Unadjusted Cox models, LumA reference. A replicated system needs concordant Basal direction, not just a small p-value.</p></div>
    <div class="figure wide">{figs['forest']}</div>
  </section>

  <section>
    <div class="section-head"><span class="section-no">06</span><h2>Survival by subtype</h2><p class="section-note">Kaplan–Meier curves. TCGA DFS and METABRIC RFS are different endpoints and are never pooled into a single survival curve.</p></div>
    <div class="figure wide">{figs['km']}</div>
  </section>

  <section>
    <div class="section-head"><span class="section-no">07</span><h2>External transport</h2><p class="section-note">The TCGA clinical + PAM50 model is locked and applied to METABRIC. Discrimination uses risk ranking; calibration at 60 months exposes baseline-hazard differences.</p></div>
    <div class="grid-2">
      <div class="figure">{figs['calibration']}</div>
      <div class="panel"><h3>What the transport test asks</h3><p>A model can rank patients well and still be miscalibrated when the baseline hazard differs. Replica reports both instead of collapsing them into one number.</p><p class="small"><strong>Discrimination.</strong> Clinical {c_base:.3f} → clinical + PAM50 {c_full:.3f}; ΔC {delta:+.3f} [{delta_lo:+.4f}, {delta_hi:+.4f}]. The interval only just excludes zero.</p><p class="small">Locked model: age per 10 years, stage III/IV, ER+, PR+, and subtype. No METABRIC outcome information enters the coefficients.</p><p class="small"><strong>Calibration slope: {transport['calibration_slope']:.2f}.</strong> A value near 1 means the locked linear predictor transfers on the log-hazard scale; a value below 1 indicates over-dispersion.</p><p class="small">{intclust_text}</p></div>
    </div>
  </section>

  <section>
    <div class="section-head"><span class="section-no">08</span><h2>Nonparametric combination</h2><p class="section-note">The four component tests are dependent. NPC calibrates the combined statistic against the joint permutation distribution instead of assuming independence.</p></div>
    <div class="grid-2">
      <div class="figure">{figs['npc']}</div>
      <div class="panel"><h3>Why not Bonferroni?</h3><p>Bonferroni controls the family-wise error rate but is conservative when tests are positively dependent. NPC preserves the empirical dependence and can detect a consistent moderate signal across OS and DFS/RFS.</p><p class="small">The global p is an association test, not a replication test. It is driven by the two METABRIC components (both at the 0.0005 floor); the TCGA components are 0.307 and 0.495. A significant global test does not rescue the cohort-specific decision.</p><p class="small">Permutation scheme: bivariate outcome pairs are permuted jointly within cohort, so censoring and the OS–DFS/RFS correlation are preserved. The global null uses {results['meta']['n_permutations_npc']:,} permutations.</p></div>
    </div>
  </section>

  <section>
    <div class="section-head"><span class="section-no">09</span><h2>Expression modules</h2><p class="section-note">NMF is fitted in TCGA only. METABRIC module scores use the locked top-30 gene lists, so external outcome information cannot leak into module definition.</p></div>
    <div class="grid-2">
      <div class="figure">{figs['rank']}</div>
      <div class="figure">{figs['modules']}</div>
    </div>
    <div class="table-wrap"><table><thead><tr><th>Module</th><th>Cohort</th><th>HR per SD</th><th>Cox p</th><th>BH q</th></tr></thead><tbody>{module_table(modules)}</tbody></table></div>
  </section>

  {simulation_section}

  <section>
    <div class="section-head"><span class="section-no">{limits_no}</span><h2>What would change the answer</h2></div>
    <div class="grid-2">
      <div class="panel"><h3>Threats to validity</h3><p>PAM50 calls are not assay-identical across cohorts. TCGA labels come from a published classifier; METABRIC labels come from a different pipeline. TCGA grade is absent from the cBioPortal export, so the primary clinical baseline uses age, stage, ER, and PR; grade is a METABRIC-only sensitivity. DFS and RFS are not the same endpoint. Complete-case analysis can bias estimates if missingness is informative.</p></div>
      <div class="panel"><h3>What the audit does not claim</h3><p>Replica is not a causal analysis, not a new classifier, and not a clinical decision tool. A subtype coefficient is an association under the fitted model. The permutation test addresses exchangeability, not confounding. The dashboard is a communication layer over frozen artifacts, not a live patient-facing prognostic device.</p></div>
    </div>
  </section>

  <section>
    <div class="section-head"><span class="section-no">{repro_no}</span><h2>Reproduce it</h2></div>
    <div class="codeblock"><code>git clone https://github.com/imSheheryar/replica<br>cd replica<br>make setup<br>make pipeline</code></div>
    <p class="small">Every artifact records its inputs, checksums, seed, and timestamp. The API exposes only frozen JSON; no patient-level data are served. Tests cover NPC calibration, FDR, harmonisation, Cox fitting, and permutation behaviour.</p>
  </section>
</main>

<footer><div class="wrap footer-grid"><div><span class="mark">Replica</span><p class="small">Built by <a href="https://github.com/imSheheryar">imSheheryar</a>. Code MIT; data governed by the source cohorts.</p></div><div class="footer-links"><a href="https://github.com/imSheheryar/replica">GitHub</a><a href="https://github.com/imSheheryar/replica/blob/main/docs/methods.md">Methods</a><a href="https://www.cbioportal.org">cBioPortal</a></div></div></footer>
</body>
</html>"""
    (SITE / "index.html").write_text(html)
    return SITE / "index.html"


if __name__ == "__main__":
    print(build())
