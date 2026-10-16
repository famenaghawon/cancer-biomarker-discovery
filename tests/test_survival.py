"""Offline tests for the Cox PH and Kaplan–Meier modules."""
import numpy as np
import pandas as pd

from cancer_biomarkers.survival import (
    cox_ph,
    kaplan_meier,
    stratify_by_expression,
)


def _toy_cohort(n=100, seed=0):
    rng = np.random.default_rng(seed)
    expr = pd.DataFrame({
        "GENE_A": rng.normal(5, 1, n),
        "GENE_B": rng.normal(5, 1, n),
        "GENE_C": rng.normal(5, 1, n),
    })
    # Time depends weakly on GENE_A only
    risk = np.exp(0.3 * (expr["GENE_A"] - 5))
    time = rng.exponential(scale=100.0 / risk)
    event = rng.binomial(1, 0.6, n)
    return expr, pd.Series(time, name="time"), pd.Series(event, name="event")


def test_cox_ph_returns_table_and_concordance():
    expr, time, event = _toy_cohort()
    res = cox_ph(expr, time, event, top_genes=["GENE_A", "GENE_B"])
    assert res.coefficients.shape[0] >= 2
    assert np.isfinite(res.concordance)
    assert 0.0 <= res.concordance <= 1.0


def test_cox_ph_rejects_zero_events():
    import pytest
    expr, time, event = _toy_cohort()
    event[:] = 0
    with pytest.raises(RuntimeError, match="Too few events"):
        cox_ph(expr, time, event, top_genes=["GENE_A", "GENE_B"])


def test_kaplan_meier_returns_monotone_curves():
    _, time, event = _toy_cohort()
    curves = kaplan_meier(time, event)
    assert len(curves) == 1
    surv = curves[0].survival
    assert np.all(np.diff(surv) <= 1e-12)   # non-increasing
    assert 0.0 <= surv.min() <= 1.0
    assert 0.0 <= surv.max() <= 1.0


def test_stratify_by_expression_median_split():
    expr = pd.Series([1, 2, 3, 4, 5, 6])
    s = stratify_by_expression(expr, cutoff="median")
    assert (s == "low").sum() == 3
    assert (s == "high").sum() == 3
