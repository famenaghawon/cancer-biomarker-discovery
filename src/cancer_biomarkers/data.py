"""Backwards-compatible shim. Prefer ``cancer_biomarkers.data_source``."""
from __future__ import annotations

from .data_source.synthetic import SyntheticDataset, generate_dataset

__all__ = ["SyntheticDataset", "generate_dataset"]
