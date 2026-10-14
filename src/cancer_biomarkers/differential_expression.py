"""Differential expression: Welch's t-test + Benjamini-Hochberg FDR."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .config import DEConfig
from .logging_utils import get_logger

logger = get_logger(__name__)

_EPS = 1e-9


def differential_expression(
    expression: pd.DataFrame,
    labels: pd.Series,
    cfg: DEConfig | None = None,
) -> pd.DataFrame:
    """Compute log2FC, Welch t-test p-values and BH-FDR adjusted q-values."""
    cfg = cfg or DEConfig()
    tumor = expression.loc[labels == 1]
    normal = expression.loc[labels == 0]

    tumor_mean = tumor.mean(axis=0)
    normal_mean = normal.mean(axis=0)

    log2fc = np.log2((tumor_mean + _EPS) / (normal_mean + _EPS))

    t_stat, p_val = stats.ttest_ind(
        tumor.values, normal.values, axis=0, equal_var=False, nan_policy="omit"
    )

    p_val = np.where(np.isnan(p_val), 1.0, p_val)
    q_val = multipletests(p_val, method=cfg.method)[1]

    df = pd.DataFrame({
        "gene": expression.columns,
        "tumor_mean": tumor_mean.values,
        "normal_mean": normal_mean.values,
        "log2fc": log2fc.values,
        "t_stat": t_stat,
        "p_value": p_val,
        "q_value": q_val,
    })

    df["significant"] = (
        (df["q_value"] < cfg.alpha) & (df["log2fc"].abs() >= cfg.min_abs_log2fc)
    )
    df["direction"] = np.where(df["log2fc"] > 0, "up", "down")
    df = df.sort_values("q_value").reset_index(drop=True)

    logger.info(
        "DE done: %d/%d significant (q<%.3f, |log2FC|>=%.2f) | "
        "top |log2FC|=%.2f, min q=%.2e",
        int(df["significant"].sum()), len(df), cfg.alpha, cfg.min_abs_log2fc,
        float(df["log2fc"].abs().max()),
        float(df["q_value"].min()),
    )
    return df