"""Centralised configuration for the pipeline."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class DataConfig:
    """Configuration for data loading and generation.

    Field defaults are tuned so the synthetic cohort is realistic but not
    trivially separable. See ``docs/methods.md`` for a discussion of
    ``effect_size`` and expected classifier performance.
    """

    source: str = "synthetic"          # "synthetic" | "tcga_brca" | "methylation"
    batch_correct: bool = False
    n_samples: int = 200
    n_genes: int = 500
    tumor_ratio: float = 0.5
    n_true_biomarkers: int = 25
    effect_size: float = 0.5
    noise_sd: float = 1.0
    batch_effect_sd: float = 0.25
    random_state: int = 42
    gene_panel: tuple[str, ...] = ()


@dataclass(frozen=True)
class DEConfig:
    """Differential-expression thresholds and multiple-testing settings."""

    alpha: float = 0.05
    min_abs_log2fc: float = 0.6
    method: str = "fdr_bh"


@dataclass(frozen=True)
class MLConfig:
    """Machine-learning split, CV, and feature-selection settings."""

    test_size: float = 0.2
    n_splits: int = 5
    top_k_features: int = 30
    random_state: int = 42


@dataclass(frozen=True)
class PathConfig:
    """Where the pipeline writes its results."""

    results_dir: Path = field(default_factory=lambda: Path("results"))
    figures_dir: Path = field(default_factory=lambda: Path("results") / "figures")

    def ensure(self) -> None:
        """Create the results and figures directories if missing."""
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)

PATHS = PathConfig()
