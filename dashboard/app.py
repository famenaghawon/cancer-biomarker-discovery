"""Interactive dashboard for the cancer-biomarkers pipeline.

Reads the artefacts under ``results/`` produced by
``python scripts/run_pipeline.py`` and renders them in the browser.

Run with:
    streamlit run dashboard/app.py
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

# ---------------------------------------------------------------------------
# Page setup
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="Cancer Biomarker Discovery",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

RESULTS = Path("results")
FIGURES = RESULTS / "figures"

# ---------------------------------------------------------------------------
# Sidebar — pick which artefacts to view
# ---------------------------------------------------------------------------

st.sidebar.title("🧬 Cancer Biomarker Discovery")
st.sidebar.markdown(
    "Interactive view over the outputs of the most recent pipeline run. "
    "Run `python scripts/run_pipeline.py` first to (re)generate them."
)

if not RESULTS.exists():
    st.error(
        f"No `{RESULTS}/` directory found. "
        "Run `python scripts/run_pipeline.py` from the repository root, then reload."
    )
    st.stop()

# ---------------------------------------------------------------------------
# Summary metrics
# ---------------------------------------------------------------------------

summary_path = RESULTS / "analysis_summary.json"
if not summary_path.exists():
    st.warning(
        "`results/analysis_summary.json` not found — showing what exists.",
        icon="⚠️",
    )
    summary = {}
else:
    summary = json.loads(summary_path.read_text())

st.title("Cancer Biomarker Discovery Pipeline")
st.caption("Synthetic and real TCGA-BRCA data · differential expression · leakage-safe ML · survival")

# Header metrics
col1, col2, col3, col4 = st.columns(4)

with col1:
    st.metric(
        "Data source",
        summary.get("config", {}).get("data", {}).get("source", "?"),
    )

with col2:
    st.metric(
        "Significant biomarkers",
        summary.get("n_significant_biomarkers", 0),
    )

with col3:
    rec = summary.get("top_k_recovery", {})
    st.metric(
        "Top-K recovery",
        f"{rec.get('n_true_in_top_k', 0)}/{rec.get('k', 0)}",
        help="How many of the top-K selected features are known ground-truth biomarkers.",
    )

with col4:
    models = summary.get("ml", {}).get("models", {})
    if models:
        best_name = max(models, key=lambda n: models[n].get("roc_auc") or 0.0)
        best_auc = models[best_name].get("roc_auc")
        st.metric(
            "Best ROC-AUC",
            f"{best_auc:.3f}" if isinstance(best_auc, (int, float)) else "n/a",
            help=f"Best model: {best_name}",
        )
    else:
        st.metric("Best ROC-AUC", "n/a")

# Survival metric (only if present)
survival = summary.get("survival")
if survival:
    c1, c2 = st.columns(2)
    with c1:
        conc = survival.get("concordance")
        st.metric(
            "Cox PH concordance",
            f"{conc:.3f}" if isinstance(conc, (int, float)) else "n/a",
        )
    with c2:
        st.metric("Survival panel", ", ".join(survival.get("panel", [])))

st.divider()

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

st.subheader("Figures")

fig_specs = [
    ("volcano.png", "Volcano plot", "log₂FC vs. −log₁₀(q), tumour vs. normal"),
    ("heatmap.png", "Heatmap", "Top biomarkers, log₂-transformed and z-scored"),
    ("roc_curve.png", "ROC curve", "Best classifier on the held-out test set"),
    ("kaplan_meier.png", "Kaplan–Meier", "Survival by high/low expression (if available)"),
]
present = [(f, t, c) for f, t, c in fig_specs if (FIGURES / f).exists()]

if not present:
    st.info("No figures found in `results/figures/`. Run the pipeline first.")
else:
    cols = st.columns(2)
    for i, (fname, title, caption) in enumerate(present):
        with cols[i % 2]:
            st.markdown(f"**{title}**")
            st.image(str(FIGURES / fname), caption=caption, width="stretch")

st.divider()

# ---------------------------------------------------------------------------
# Differential expression table
# ---------------------------------------------------------------------------

st.subheader("Differential expression")

de_path = RESULTS / "de_results.csv"
if de_path.exists():
    de = pd.read_csv(de_path)

    # Quick stats
    c1, c2, c3 = st.columns(3)
    with c1:
        st.metric("Total features tested", len(de))
    with c2:
        st.metric("Significant (q < 0.05, |log₂FC| ≥ 0.6)", int(de["significant"].sum()))
    with c3:
        st.metric("Max |log₂FC|", f"{de['log2fc'].abs().max():.2f}")

    # Interactive volcano
    vol = de.copy()
    vol["-log10q"] = -vol["q_value"].clip(lower=1e-300).apply(lambda x: -__import__("math").log10(x))
    # Simpler: use numpy
    import numpy as np
    vol["-log10q"] = -np.log10(vol["q_value"].clip(lower=1e-300))
    fig = px.scatter(
        vol,
        x="log2fc",
        y="-log10q",
        color="significant",
        hover_data=["gene", "p_value", "q_value"],
        title="Volcano (interactive)",
        labels={"log2fc": "log₂ fold change", "-log10q": "−log₁₀(q)"},
    )
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="grey")
    fig.add_vline(x=0.6, line_dash="dash", line_color="grey")
    fig.add_vline(x=-0.6, line_dash="dash", line_color="grey")
    st.plotly_chart(fig, width="stretch")

    # Table with sortable columns
    st.markdown("**Top 50 by q-value**")
    st.dataframe(
        de.sort_values("q_value").head(50),
        width="stretch",
        height=400,
    )

    # Download button
    st.download_button(
        "Download full DE table (CSV)",
        data=de.to_csv(index=False).encode(),
        file_name="de_results.csv",
        mime="text/csv",
    )
else:
    st.info("No `de_results.csv` found. Run the pipeline first.")

st.divider()

# ---------------------------------------------------------------------------
# ML model metrics
# ---------------------------------------------------------------------------

st.subheader("Machine learning")

if models:
    rows = []
    for name, m in models.items():
        rows.append({
            "model": name,
            "accuracy": round(m.get("accuracy", 0.0), 3),
            "precision": round(m.get("precision", 0.0), 3),
            "recall": round(m.get("recall", 0.0), 3),
            "f1": round(m.get("f1", 0.0), 3),
            "roc_auc": round(m.get("roc_auc") or 0.0, 3),
            "cv_roc_auc_mean": round(m.get("cv_roc_auc_mean", 0.0), 3),
            "cv_roc_auc_std": round(m.get("cv_roc_auc_std", 0.0), 3),
        })
    metrics_df = pd.DataFrame(rows).set_index("model")
    st.dataframe(metrics_df, width="stretch")

    # Bar chart
    fig = px.bar(
        metrics_df.reset_index(),
        x="model",
        y=["accuracy", "f1", "roc_auc"],
        barmode="group",
        title="Model comparison",
    )
    st.plotly_chart(fig, width="stretch")
else:
    st.info("No `ml` block in the summary JSON.")

# ---------------------------------------------------------------------------
# Survival (if present)
# ---------------------------------------------------------------------------

if survival:
    st.divider()
    st.subheader("Survival")
    st.markdown(
        f"**Cox PH concordance:** {survival.get('concordance', 'n/a')}  \n"
        f"**Panel:** {', '.join(survival.get('panel', []))}  \n"
        f"**Top gene for KM:** {survival.get('top_gene', 'n/a')}"
    )
    cox_path = RESULTS / "cox_coefficients.csv"
    if cox_path.exists():
        cox = pd.read_csv(cox_path, index_col=0)
        st.markdown("**Cox PH coefficients**")
        st.dataframe(cox, width="stretch")

# ---------------------------------------------------------------------------
# Footer
# ---------------------------------------------------------------------------

st.divider()
st.caption(
    "Synthetic and public TCGA-BRCA data only. "
    "Not for clinical use. See `docs/methods.md` for full methodology."
)
