"""Data sources: synthetic expression, TCGA-BRCA, and synthetic methylation."""
from .synthetic import SyntheticDataset, generate_dataset
from .tcga import TCGADataset, load_tcga_brca

__all__ = [
    "SyntheticDataset",
    "TCGADataset",
    "generate_dataset",
    "load_tcga_brca",
]
