"""Publication-style figures."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .logging_utils import get_logger

logger = get_logger(__name__)
sns.set_theme(style="whitegrid")


def volcano_plot(de_df: pd.DataFrame, out: Path) -> None:
    df = de_df.copy()
    df["-log10q"] = -np.log10(df["q_value"].clip(lower=1e-300))
    sig = df["significant"]

    plt.figure(figsize=(8, 6))
    plt.scatter(df.loc[~sig, "log2fc"], df.loc[~sig, "-log10q"],
                s=8, alpha=0.4, label="n.s.", color="grey")
    plt.scatter(df.loc[sig & (df["log2fc"] > 0), "log2fc"],
                df.loc[sig & (df["log2fc"] > 0), "-log10q"],
                s=14, alpha=0.8, label="up", color="crimson")
    plt.scatter(df.loc[sig & (df["log2fc"] < 0), "log2fc"],
                df.loc[sig & (df["log2fc"] < 0), "-log10q"],
                s=14, alpha=0.8, label="down", color="steelblue")
    plt.axhline(-np.log10(0.05), ls="--", c="k", lw=0.8)
    plt.axvline(1, ls="--", c="k", lw=0.8); plt.axvline(-1, ls="--", c="k", lw=0.8)
    plt.xlabel("log2 fold change"); plt.ylabel("-log10(q)")
    plt.title("Volcano plot — tumor vs normal")
    plt.legend(); plt.tight_layout()
    plt.savefig(out, dpi=150); plt.close()
    logger.info("Saved %s", out)


def heatmap_top_genes(expression: pd.DataFrame, labels: pd.Series,
                      genes: list[str], out: Path) -> None:
    sub = expression[genes].copy()
    sub["label"] = labels.values
    sub = sub.sort_values("label")
    mat = sub.drop(columns="label").T
    plt.figure(figsize=(10, max(4, 0.2 * len(genes))))
    sns.heatmap(mat, cmap="vlag", center=mat.values.mean(),
                xticklabels=False, yticklabels=True, cbar_kws={"label": "expr"})
    plt.title("Top biomarkers — tumor vs normal")
    plt.tight_layout(); plt.savefig(out, dpi=150); plt.close()
    logger.info("Saved %s", out)
