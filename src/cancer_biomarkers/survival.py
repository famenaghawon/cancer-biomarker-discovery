"""Survival analysis on the loaded cohort.

Cox proportional-hazards fit is provided by ``lifelines.CoxPHFitter``;
Kaplan–Meier curves are computed via a small product-limit estimator that
groups deaths by unique time (see ``kaplan_meier`` below).

Covariate hygiene
-----------------
Cox PH fails when covariates are collinear or degenerate. Because our top-K
feature selection often picks highly co-expressed genes (e.g. several
proliferation markers), we drop near-duplicate covariates before fitting.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class CoxResult:
    coefficients: pd.DataFrame       # gene -> coef, HR, CI, p-value
    concordance: float


@dataclass
class KMCurve:
    times: np.ndarray
    survival: np.ndarray
    label: str


_CORR_THRESHOLD = 0.98


def _filter_collinear(X: pd.DataFrame, cols: list[str]) -> list[str]:
    """Greedily drop covariates correlated above ``_CORR_THRESHOLD``."""
    keep = [c for c in cols if X[c].var() > 1e-8]
    if len(keep) < 2:
        return keep
    kept: list[str] = []
    for col in keep:
        if not kept:
            kept.append(col)
            continue
        corr = X[kept + [col]].corr().abs().loc[col, kept].max()
        if corr < _CORR_THRESHOLD:
            kept.append(col)
    return kept


def cox_ph(
    expression: pd.DataFrame,
    time: pd.Series,
    event: pd.Series,
    top_genes: list[str],
) -> CoxResult:
    """Fit a Cox PH model via lifelines.

    Rows with missing time, event, or gene expression are dropped.
    Covariates are standardized so hazard ratios are per 1 SD of expression.
    """
    from lifelines import CoxPHFitter

    panel = [g for g in top_genes if g in expression.columns]
    if len(panel) < 2:
        raise ValueError(f"Need at least 2 covariates, got {panel}")

    df = pd.concat(
        [expression[panel], time.rename("time"), event.rename("event")],
        axis=1,
    ).dropna(subset=["time", "event", *panel])
    df = df[df["time"] > 0]

    if len(df) < 20:
        raise RuntimeError(f"Too few usable records ({len(df)})")
    n_events = int(df["event"].sum())
    if n_events < 10:
        raise RuntimeError(f"Too few events ({n_events})")

    # Standardize covariates; drop zero-variance columns.
    raw = df[panel].astype(float)
    keep = [c for c in panel if raw[c].var() > 1e-8]
    if len(keep) < 2:
        raise RuntimeError(f"Fewer than 2 non-degenerate covariates: {keep}")
    std = raw[keep].std().replace(0, 1.0)
    df_std = (raw[keep] - raw[keep].mean()) / std

    fit_df = pd.concat(
        [df_std, df[["time", "event"]]], axis=1,
    ).rename(columns={"time": "T", "event": "E"})

    cph = CoxPHFitter()
    cph.fit(fit_df, duration_col="T", event_col="E")

    summary = cph.summary
    table = pd.DataFrame({
        "coef": summary["coef"],
        "hazard_ratio": summary["exp(coef)"],
        "hr_ci_low": summary["exp(coef) lower 95%"],
        "hr_ci_high": summary["exp(coef) upper 95%"],
        "p_value": summary["p"],
    })
    table.attrs["interpretation"] = "HR per 1 SD of standardized gene expression"

    concordance = float(cph.concordance_index_)
    logger.info(
        "Cox PH: %d samples, %d events, %d covariates, concordance=%.3f",
        len(fit_df), n_events, len(keep), concordance,
    )
    return CoxResult(coefficients=table, concordance=concordance)


def kaplan_meier(
    time: pd.Series,
    event: pd.Series,
    stratify_by: pd.Series | None = None,
    groups: tuple[str, str] = ("low", "high"),
) -> list[KMCurve]:
    """Compute Kaplan-Meier survival curves via ``lifelines.KaplanMeierFitter``.

    Tied event times are handled by the standard product-limit estimator.
    """
    from lifelines import KaplanMeierFitter

    df = pd.DataFrame({"time": time, "event": event})
    if stratify_by is not None:
        df["stratum"] = stratify_by

    strata = (
        [(g, df[df["stratum"] == g]) for g in groups]
        if stratify_by is not None
        else [("all", df)]
    )
    curves: list[KMCurve] = []
    for label, sub in strata:
        sub = sub.dropna(subset=["time", "event"])
        sub = sub[sub["time"] > 0]
        if sub.empty:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(sub["time"], event_observed=sub["event"], label=label)
        curves.append(KMCurve(
            times=kmf.survival_function_.index.to_numpy(dtype=float),
            survival=kmf.survival_function_[label].to_numpy(dtype=float),
            label=label,
        ))
    return curves

def stratify_by_expression(
    expression: pd.Series, cutoff: str = "median"
) -> pd.Series:
    """Split a gene's expression into 'low' / 'high' at the given cutoff."""
    if cutoff == "median":
        thr = expression.median()
    elif cutoff == "tertile_low":
        thr = expression.quantile(1 / 3)
    elif cutoff == "tertile_high":
        thr = expression.quantile(2 / 3)
    else:
        raise ValueError(f"Unknown cutoff: {cutoff!r}")
    return pd.Series(
        np.where(expression >= thr, "high", "low"),
        index=expression.index,
        name="stratum",
    )
