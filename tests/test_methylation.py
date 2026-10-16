"""Methylation generator and M-value transform."""
import pandas as pd

from cancer_biomarkers.config import DataConfig
from cancer_biomarkers.multiomics import generate_methylation, to_m_values


def test_methylation_beta_range():
    cfg = DataConfig(n_samples=60, n_genes=80, n_true_biomarkers=8)
    ds = generate_methylation(cfg)
    assert ds.beta.shape == (60, 80)
    assert ds.beta.min().min() >= 0.0
    assert ds.beta.max().max() <= 1.0
    assert len(ds.true_biomarkers) == 8


def test_methylation_uses_real_cpg_ids():
    cfg = DataConfig(n_samples=30, n_genes=50, n_true_biomarkers=5)
    ds = generate_methylation(cfg)
    # Every probe should look like a real 450k ID: cg########
    assert all(col.startswith("cg") and len(col) == 10 for col in ds.beta.columns)


def test_m_value_transform_inverts_beta_ordering():
    # β = 0.5 → M = 0; β > 0.5 → M > 0; β < 0.5 → M < 0
    beta = pd.DataFrame({"cg00000001": [0.1, 0.5, 0.9]})
    m = to_m_values(beta)
    assert m.iloc[0, 0] < 0
    assert abs(m.iloc[1, 0]) < 1e-6
    assert m.iloc[2, 0] > 0
