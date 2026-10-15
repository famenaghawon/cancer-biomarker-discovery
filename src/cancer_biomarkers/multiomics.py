"""Multi-omics stub: methylation beta values.

Not integrated with the main pipeline yet. Provides an interface that
matches ``data_source.synthetic`` so a future extension can plug in a
methylation matrix and reuse the DE + ML stack.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .config import DataConfig
from .logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class MethylationDataset:
    beta: pd.DataFrame           # samples x CpG sites, values in [0, 1]
    labels: pd.Series
    clinical: pd.DataFrame
    true_biomarkers: list[str]


def generate_methylation(cfg: DataConfig) -> MethylationDataset:
    """Synthetic methylation beta values with a tumor-vs-normal signal.

    Beta values follow a bimodal distribution (unmethylated ~0.05,
    methylated ~0.95), which is what real 450k arrays look like.
    """
    rng = np.random.default_rng(cfg.random_state)
    n, g = cfg.n_samples, cfg.n_genes
    n_tumor = int(round(n * cfg.tumor_ratio))  # noqa: RUF046
    labels = np.array([1] * n_tumor + [0] * (n - n_tumor))
    rng.shuffle(labels)

    # Baseline: each site is unmethylated in normal tissue
    beta = rng.beta(0.5, 9.5, size=(n, g)).astype(float)

    # Tumor-specific hypermethylation at a subset of sites
    n_true = cfg.n_true_biomarkers
    site_names = [f"cg{i:05d}" for i in range(g)]
    true_sites = rng.choice(site_names, size=n_true, replace=False)
    true_idx = [site_names.index(s) for s in true_sites]
    for i in true_idx:
        shift = rng.uniform(0.3, 0.6)
        beta[labels == 1, i] = np.clip(beta[labels == 1, i] + shift, 0, 1)

    beta_df = pd.DataFrame(beta, columns=site_names)
    labels_s = pd.Series(labels, name="label")
    clinical = pd.DataFrame({
        "sample_id": [f"S{i:04d}" for i in range(n)],
        "label": labels,
    })
    logger.info("Generated methylation matrix: %d x %d", n, g)
    return MethylationDataset(beta_df, labels_s, clinical, list(true_sites))
