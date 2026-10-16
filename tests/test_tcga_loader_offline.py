"""Offline tests for the TCGA clinical loader.

These tests do not hit the network. They exercise ``_build_labels`` and
``_build_clinical_summary`` on synthetic frames that mimic the Xena schema.
"""
import numpy as np
import pandas as pd

from cancer_biomarkers.data_source.tcga import (
    _build_clinical_summary,
    _build_labels,
)


def _fake_xena_clinical(n_samples=20):
    idx = [f"S{i:04d}" for i in range(n_samples)]
    sample_types = (["Primary Tumor"] * (n_samples - 4)
                    + ["Solid Tissue Normal"] * 4)
    vital = ["DECEASED"] * 6 + ["LIVING"] * (n_samples - 6)
    dtd = [100.0] * 6 + [np.nan] * (n_samples - 6)
    dtf = [np.nan] * 6 + [500.0] * (n_samples - 6)
    return pd.DataFrame({
        "sample_type": sample_types,
        "pathologic_stage": ["Stage IIA"] * n_samples,
        "vital_status": vital,
        "days_to_death": dtd,
        "days_to_last_followup": dtf,
        "age_at_initial_pathologic_diagnosis": [55] * n_samples,
        "gender": ["FEMALE"] * n_samples,
    }, index=idx)


def test_build_labels_separates_tumour_and_normal():
    clinical = _fake_xena_clinical()
    labels = _build_labels(clinical)
    assert labels.sum() == 16    # 16 primary tumors
    assert (labels == 0).sum() == 4   # 4 solid tissue normals


def test_build_clinical_summary_creates_event_and_time():
    clinical = _fake_xena_clinical()
    labels = _build_labels(clinical)
    summary = _build_clinical_summary(clinical, labels)
    assert "event" in summary.columns
    assert "time_days" in summary.columns
    assert summary["event"].sum() == 6
    # Deceased rows use days_to_death, others use days_to_last_followup
    deceased = summary.loc[summary["event"] == 1]
    assert (deceased["time_days"] == 100.0).all()
    censored = summary.loc[summary["event"] == 0]
    assert (censored["time_days"] == 500.0).all()
