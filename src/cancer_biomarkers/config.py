"""Centralised configuration for the pipeline."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class DataConfig:
    n_samples: int = 200
    n_genes: int = 500
    tumor_ratio: float = 0.5
    n_true_biomarkers: int = 25
    effect_size: float = 0.5          # log2-scale mean shift for true biomarkers
    noise_sd: float = 0.75            # kept for backward compat
    batch_effect_sd: float = 0.25     # kept for backward compat
    random_state: int = 42


@dataclass(frozen=True)
class DEConfig:
    alpha: float = 0.05
    min_abs_log2fc: float = 0.8
    method: str = "fdr_bh"


@dataclass(frozen=True)
class MLConfig:
    test_size: float = 0.2
    n_splits: int = 5
    top_k_features: int = 30
    random_state: int = 42


@dataclass(frozen=True)
class PathConfig:
    results_dir: Path = field(default_factory=lambda: Path("results"))
    figures_dir: Path = field(default_factory=lambda: Path("results") / "figures")

    def ensure(self) -> None:
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)


PATHS = PathConfig()
