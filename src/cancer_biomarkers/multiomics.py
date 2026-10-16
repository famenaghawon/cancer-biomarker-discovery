"""Multi-omics: synthetic methylation data with real CpG probe labels.

The generator produces β-values on [0, 1] — the standard output of Illumina
450k arrays. β is bimodal (unmethylated ~0.05, methylated ~0.95) in real
tissue; we reproduce that, and add a tumor-specific hypermethylation signal
on a subset of the curated CpG panel.

The M-value transform in ``to_m_values`` converts β to a log-odds scale so
the downstream DE + ML pipeline (which expects log-scale continuous values)
can be reused unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .config import DataConfig
from .data_source.cpg_panel import flat_probe_list
from .logging_utils import get_logger

logger = get_logger(__name__)

_EPS = 1e-6


@dataclass
class MethylationDataset:
    beta: pd.DataFrame           # samples x CpG sites, values in [0, 1]
    labels: pd.Series            # 0 = normal, 1 = tumor
    clinical: pd.DataFrame
    true_biomarkers: list[str]   # CpG IDs that carry a real tumor-vs-normal signal


def _build_probe_names(n_probes: int, rng: np.random.Generator) -> list[str]:
    """Prefer real curated CpG IDs; fill remainder with plausible fake IDs."""
    curated = flat_probe_list()
    if n_probes <= len(curated):
        picks = list(rng.choice(curated, size=n_probes, replace=False))
        rng.shuffle(picks)
        return picks

    # Fill after curated pool is exhausted (for large n_genes configs).
    picks = list(curated)
    remaining = n_probes - len(curated)
    fake = [f"cg{90000000 + i:08d}" for i in range(remaining)]
    all_probes = picks + fake
    rng.shuffle(all_probes)
    return all_probes


def generate_methylation(cfg: DataConfig) -> MethylationDataset:
    """Generate a synthetic methylation β-matrix with a tumor signal."""
    rng = np.random.default_rng(cfg.random_state)
    n, g = cfg.n_samples, cfg.n_genes
    n_tumor = int(round(n * cfg.tumor_ratio))  # noqa: RUF046
    labels = np.array([1] * n_tumor + [0] * (n - n_tumor))
    rng.shuffle(labels)

    # Baseline β: mostly unmethylated, with a tail toward hypermethylation
    beta = rng.beta(a=0.5, b=9.5, size=(n, g)).astype(float)

    probe_names = _build_probe_names(g, rng)
    n_true = min(cfg.n_true_biomarkers, g)
    true_probes = rng.choice(probe_names, size=n_true, replace=False).tolist()
    probe_index = {p: i for i, p in enumerate(probe_names)}

    # 15% of tumours behave like normals on the true panel (matching the
    # expression-generator design so classifier performance stays realistic).
    silent = (labels == 1) & (rng.random(n) < 0.15)
    tumor_mask = (labels == 1) & (~silent)

    for probe in true_probes:
        i = probe_index[probe]
        shift = rng.uniform(0.05, 0.15)   # moderate hypermethylation shift
        beta[tumor_mask, i] = np.clip(beta[tumor_mask, i] + shift, 0.0, 1.0)
        beta[silent, i] = np.clip(
            beta[silent, i] + rng.normal(0, 0.02, size=int(silent.sum())),
            0.0, 1.0,
        )

    # Global noise floor; keep beta strictly inside (0, 1) for the M-value
    # transform to be well-defined.
    beta = beta + rng.normal(0, 0.02, size=beta.shape)
    beta = np.clip(beta, 1e-3, 1.0 - 1e-3)

    beta_df = pd.DataFrame(beta, columns=probe_names)
    labels_s = pd.Series(labels, name="label")
    clinical = pd.DataFrame({
        "sample_id": [f"S{i:04d}" for i in range(n)],
        "label": labels,
        "batch": rng.choice([0, 1, 2], size=n),
    })

    logger.info(
        "Generated methylation: %d samples x %d CpG sites (%d with tumor signal, "
        "%d silent tumours, %d from curated 450k panel)",
        n, g, len(true_probes), int(silent.sum()),
        min(g, len(flat_probe_list())),
    )
    return MethylationDataset(beta_df, labels_s, clinical, true_probes)


def to_m_values(beta: pd.DataFrame) -> pd.DataFrame:
    """Convert β-values to M-values (log2-odds), matching the DE input scale.

    ``M = log2((β + ε) / (1 - β + ε))`` with ε = 1e-6 to avoid infinities at
    the boundaries. This is the standard transformation for array-based
    methylation analysis; see Du et al. (2010), BMC Bioinformatics 11:587.
    """
    # Clip in numpy to avoid any pandas clipping edge cases on 2-D frames.
    arr = beta.to_numpy(dtype=float)
    arr = np.clip(arr, _EPS, 1.0 - _EPS)
    m = np.log2(arr / (1.0 - arr))
    return pd.DataFrame(m, index=beta.index, columns=beta.columns)
