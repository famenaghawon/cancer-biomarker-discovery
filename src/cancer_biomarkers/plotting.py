"""Publication-style figures."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import roc_curve

from .logging_utils import get_logger

logger = get_logger(__name__)
sns.set_theme(style="whitegrid")


def volcano_plot(
    de_df: pd.DataFrame,
    out: Path,
    min_abs_log2fc: float = 0.6,
    alpha: float = 0.05,
) -> None:
    """Volcano plot with threshold lines matching the DE significance rule."""
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

    plt.axhline(-np.log10(alpha), ls="--", c="k", lw=0.8)
    plt.axvline(min_abs_log2fc, ls="--", c="k", lw=0.8)
    plt.axvline(-min_abs_log2fc, ls="--", c="k", lw=0.8)

    plt.xlabel("log2 fold change")
    plt.ylabel("-log10(q)")
    plt.title(f"Volcano plot — tumor vs normal (|log2FC| >= {min_abs_log2fc}, q < {alpha})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()
    logger.info("Saved %s", out)


def heatmap_top_genes(
    expression: pd.DataFrame,
    labels: pd.Series,
    genes: list[str],
    out: Path,
    log_transform: bool = True,
) -> None:
    """Heatmap of top biomarkers, optionally log2-transformed and z-scored per gene."""
    sub = expression[genes].copy()
    if log_transform:
        sub = np.log2(sub.clip(lower=1e-9))
    # z-score each gene (row of mat) so patterns are visually comparable
    mat = sub.T
    mat = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, 1), axis=0)

    # Order samples by class then by first gene's z-score
    order = labels.sort_values(kind="stable").index
    mat = mat.loc[:, order]

    plt.figure(figsize=(10, max(4, 0.2 * len(genes))))
    sns.heatmap(
        mat, cmap="vlag", center=0.0,
        xticklabels=False, yticklabels=True,
        cbar_kws={"label": "z-score (per gene)"},
    )
    plt.title("Top biomarkers — tumor vs normal (z-scored, log2)")
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()
    logger.info("Saved %s", out)


def roc_curve_plot(
    y_true: np.ndarray,
    y_score: np.ndarray,
    out: Path,
    label: str = "model",
) -> None:
    """ROC curve for a binary classifier."""
    fpr, tpr, _ = roc_curve(y_true, y_score)
    from sklearn.metrics import roc_auc_score
    auc = roc_auc_score(y_true, y_score)

    plt.figure(figsize=(6, 6))
    plt.plot(fpr, tpr, lw=2, label=f"{label} (AUC = {auc:.3f})")
    plt.plot([0, 1], [0, 1], ls="--", c="k", lw=0.8)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title("ROC curve")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()
    logger.info("Saved %s", out)
