"""The recovery metric must be omitted for sources without a planted panel."""
import numpy as np
import pandas as pd

from cancer_biomarkers.config import DataConfig, PathConfig


class _FakeNoPanelDataset:
    """Minimal stand-in for a TCGADataset-like object with no planted panel."""

    def __init__(self, n_samples=60, n_features=40, seed=0):
        rng = np.random.default_rng(seed)
        n_tumor = n_samples // 2
        self.expression = pd.DataFrame(
            rng.normal(8, 1, size=(n_samples, n_features)),
            columns=[f"G{i:03d}" for i in range(n_features)],
        )
        # Plant a small, real signal so DE has something to find.
        labels = np.array([1] * n_tumor + [0] * (n_samples - n_tumor))
        rng.shuffle(labels)
        self.expression.loc[labels == 1, "G000"] += 1.5
        self.labels = pd.Series(labels, name="label")
        self.clinical = pd.DataFrame({"sample_id": [f"S{i:03d}" for i in range(n_samples)]})
        # No planted panel: this is what real TCGA looks like.
        self.gene_panel = list(self.expression.columns)
        self.has_planted_biomarkers = False


def test_recovery_omitted_when_no_planted_panel(tmp_path, monkeypatch):
    """A dataset with ``has_planted_biomarkers=False`` must not emit ``top_k_recovery``."""
    from cancer_biomarkers import pipeline as pipeline_module

    fake = _FakeNoPanelDataset()
    monkeypatch.setattr(pipeline_module, "_load_dataset", lambda cfg: fake)

    paths = PathConfig(results_dir=tmp_path, figures_dir=tmp_path / "figures")
    summary = pipeline_module.run(
        data_cfg=DataConfig(n_samples=60, n_genes=40, n_true_biomarkers=5),
        paths=paths,
    )

    assert "top_k_recovery" not in summary, (
        "Recovery metric must be omitted when the source has no planted panel"
    )
    assert summary["true_biomarkers"] == []
    assert summary["true_biomarkers_recovered"] == []


def test_recovery_present_when_planted_panel(tmp_path):
    """Sanity check: the synthetic source still emits ``top_k_recovery``."""
    from cancer_biomarkers.pipeline import run as pipeline_run

    paths = PathConfig(results_dir=tmp_path, figures_dir=tmp_path / "figures")
    summary = pipeline_run(
        data_cfg=DataConfig(n_samples=60, n_genes=80, n_true_biomarkers=5),
        paths=paths,
    )

    assert "top_k_recovery" in summary
    assert summary["top_k_recovery"]["k"] > 0
