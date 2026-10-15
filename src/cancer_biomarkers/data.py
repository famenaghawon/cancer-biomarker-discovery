"""Synthetic TCGA-like expression data generation."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .config import DataConfig
from .logging_utils import get_logger

logger = get_logger(__name__)


@dataclass
class SyntheticDataset:
    expression: pd.DataFrame
    labels: pd.Series
    clinical: pd.DataFrame
    true_biomarkers: list[str]


_REAL_GENE_NAMES = [
    "BRCA1", "BRCA2", "TP53", "EGFR", "MYC", "KRAS", "PTEN", "RB1",
    "PIK3CA", "AKT1", "MTOR", "CDH1", "APC", "VHL", "NF1", "SMAD4",
    "ALK", "BRAF", "RET", "MET", "ERBB2", "FGFR1", "FGFR2", "FGFR3",
    "CCND1", "CDK4", "CDK6", "MDM2", "BCL2", "BAX",
]


def _make_gene_names(n_genes: int, n_true: int, rng: np.random.Generator) -> list[str]:
    real = list(_REAL_GENE_NAMES)
    rng.shuffle(real)
    true_names = real[:n_true]
    filler = [f"Gene_{i:04d}" for i in range(1, n_genes - n_true + 1)]
    names = true_names + filler
    rng.shuffle(names)
    return names


def generate_dataset(cfg: DataConfig) -> SyntheticDataset:
    rng = np.random.default_rng(cfg.random_state)
    n = cfg.n_samples
    g = cfg.n_genes

    gene_names = _make_gene_names(g, cfg.n_true_biomarkers, rng)
    true_biomarkers = [gn for gn in gene_names if not gn.startswith("Gene_")]
    idx_true = np.array([gene_names.index(b) for b in true_biomarkers])

    # Baseline expression on log2 scale
    baseline_log2 = rng.normal(loc=8.0, scale=1.0, size=g)
    per_gene_sd = rng.uniform(0.4, 1.2, size=g)

    # Labels
    n_tumor = int(round(n * cfg.tumor_ratio))  # noqa: RUF046
    labels = np.array([1] * n_tumor + [0] * (n - n_tumor))
    rng.shuffle(labels)

    # Log2 expression matrix
    noise = rng.normal(0, cfg.noise_sd, size=(n, g)) * per_gene_sd[None, :]
    log2_expr = baseline_log2[None, :] + noise

    # Per-sample library-size and batch effects, driven by cfg.batch_effect_sd
    lib = rng.normal(0, cfg.batch_effect_sd * 0.8, size=(n, 1))
    batch = rng.choice(
        [-cfg.batch_effect_sd, 0.0, cfg.batch_effect_sd],
        size=(n, 1),
    )
    log2_expr = log2_expr + lib + batch

    # Tumor effect on true biomarkers
    raw_effects = rng.normal(cfg.effect_size, cfg.effect_size * 0.5,
                             size=len(idx_true))
    raw_effects = np.clip(raw_effects, 0.2, cfg.effect_size * 2.5)
    flip = rng.random(len(idx_true)) < 0.25
    effects = raw_effects * np.where(flip, -1.0, 1.0)

    silent = (labels == 1) & (rng.random(n) < 0.20)
    tumor_mask = (labels == 1) & (~silent)
    for i, gi in enumerate(idx_true):
        log2_expr[tumor_mask, gi] += effects[i]
        log2_expr[silent, gi] += rng.normal(0, 0.15, size=int(silent.sum()))

    # Decoy genes
    n_decoys = max(10, cfg.n_true_biomarkers)
    non_true_idx = np.array(
        [i for i in range(g) if i not in set(idx_true.tolist())]
    )
    rng.shuffle(non_true_idx)
    decoys = non_true_idx[:n_decoys]
    for gi in decoys:
        decoy_effect = rng.normal(0.0, cfg.effect_size * 0.8)
        log2_expr[labels == 1, gi] += decoy_effect

    # Global noise floor
    log2_expr += rng.normal(0, cfg.noise_sd * 0.5, size=(n, g))

    # Convert to linear expression units
    expr = np.power(2.0, log2_expr)

    expression = pd.DataFrame(expr, columns=gene_names)
    labels_s = pd.Series(labels, name="label")

    # Clinical metadata: stages only for tumors, "NA" for normals
    stages = np.array(["NA"] * n, dtype=object)
    tumor_idx = np.where(labels == 1)[0]
    stages[tumor_idx] = rng.choice(
        ["I", "II", "III", "IV"],
        size=len(tumor_idx),
        p=[0.25, 0.30, 0.30, 0.15],
    )

    clinical = pd.DataFrame({
        "sample_id": [f"S{i:04d}" for i in range(n)],
        "label": labels,
        "age": rng.integers(35, 85, size=n),
        "stage": stages,
        "batch": batch.ravel(),
    })

    logger.info(
        "Generated dataset: %d samples x %d genes (%d true biomarkers, "
        "effect log2 ~ N(%.2f, %.2f), %d silent tumors, %d decoys)",
        n, g, len(true_biomarkers),
        cfg.effect_size, cfg.effect_size * 0.5,
        int(silent.sum()), n_decoys,
    )
    return SyntheticDataset(expression, labels_s, clinical, true_biomarkers)
