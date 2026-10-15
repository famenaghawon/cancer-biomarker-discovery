"""Survival analysis on the loaded cohort.

Implements a Cox proportional-hazards fit via ``statsmodels`` and a
Kaplan-Meier curve via a small hand-rolled implementation (no lifelines
dependency). Both are honest about right-censoring.

Covariate hygiene
-----------------
Cox PH fails when covariates are collinear (singular Hessian). Because our
top-K feature selection often picks highly co-expressed genes (e.g. several
proliferation markers), we drop near-duplicate covariates before fitting.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg

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
    """Fit a Cox PH model on the top-k gene expressions.

    Rows with missing time or event are dropped. Collinear covariates
    (correlation above 0.98) are removed before fitting.
    """
    df = pd.concat(
        [expression[top_genes], time.rename("time"), event.rename("event")],
        axis=1,
    )
    df = df.dropna(subset=["time", "event"])
    df = df[df["time"] > 0]
    if len(df) < 20:
        raise RuntimeError(
            f"Not enough usable survival records ({len(df)}). "
            "Check the clinical table for time/event columns."
        )
    n_events = int(df["event"].sum())
    if n_events < 10:
        raise RuntimeError(
            f"Too few events for Cox PH ({n_events}). "
            "The clinical table may use a different vital-status encoding, "
            "or the cohort may have no recorded deaths."
        )

    raw = df[top_genes].astype(float)
    # Standardize covariates: Cox PH is scale-sensitive and unstandardized
    # expression values (order 1 to 1e3) break the numerical Hessian.
    mean = raw.mean(axis=0)
    std = raw.std(axis=0).replace(0, 1.0)
    X_std = (raw - mean) / std

    kept = _filter_collinear(X_std, top_genes)
    dropped = sorted(set(top_genes) - set(kept))
    if dropped:
        logger.info(
            "Cox PH: dropped %d collinear/degenerate covariates: %s",
            len(dropped), ", ".join(dropped),
        )
    if len(kept) < 2:
        raise RuntimeError(
            "Too few non-degenerate covariates for Cox PH after filtering "
            f"(kept {len(kept)} of {len(top_genes)})"
        )

    X_kept = sm.add_constant(X_std[kept])
    y = df["time"].astype(float)
    cens = 1 - df["event"].astype(int)
    model = PHReg(y, X_kept, status=cens, ties="efron")
    fit = model.fit(method="bfgs")

    params = fit.params
    bse = fit.bse
    pvalues = fit.pvalues
    hazard_ratios = np.exp(np.clip(params, -50, 50))
    ci_low = np.exp(np.clip(params - 1.96 * bse, -50, 50))
    ci_high = np.exp(np.clip(params + 1.96 * bse, -50, 50))

    table = pd.DataFrame({
        "coef": params,
        "hazard_ratio": hazard_ratios,
        "hr_ci_low": ci_low,
        "hr_ci_high": ci_high,
        "p_value": pvalues,
    })
    table.attrs["interpretation"] = "HR per 1 SD of gene expression"

    try:
        concordance = float(fit.concordance)
    except AttributeError:
        concordance = float("nan")

    logger.info(
        "Cox PH: %d samples, %d events, %d covariates, concordance=%.3f",
        len(df), int(df["event"].sum()), len(kept), concordance,
    )
    return CoxResult(coefficients=table, concordance=concordance)


def kaplan_meier(
    time: pd.Series,
    event: pd.Series,
    stratify_by: pd.Series | None = None,
    groups: tuple[str, str] = ("low", "high"),
) -> list[KMCurve]:
    """Compute Kaplan-Meier survival curves, optionally by a binary stratum."""
    df = pd.DataFrame({"time": time, "event": event})
    if stratify_by is not None:
        df["stratum"] = stratify_by

    curves: list[KMCurve] = []
    strata = (
        [(g, df[df["stratum"] == g]) for g in groups]
        if stratify_by is not None
        else [("all", df)]
    )
    for label, sub in strata:
        sub = sub.dropna(subset=["time", "event"])
        sub = sub[sub["time"] > 0].sort_values("time")
        if sub.empty:
            continue
        times = sub["time"].to_numpy(dtype=float)
        events = sub["event"].to_numpy(dtype=int)
        surv = []
        s = 1.0
        for t, e in zip(times, events, strict=False):
            if e == 1:
                at_risk = int((times >= t).sum())
                if at_risk > 0:
                    s *= 1.0 - 1.0 / at_risk
            surv.append(s)
        curves.append(KMCurve(times=times, survival=np.array(surv), label=label))
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
