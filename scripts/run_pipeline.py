"""CLI entry point for the cancer-biomarkers pipeline.

Usage
-----
python scripts/run_pipeline.py                    # synthetic (default)
python scripts/run_pipeline.py --data tcga_brca   # real TCGA-BRCA subset
python scripts/run_pipeline.py --nested           # also run nested CV
"""
from __future__ import annotations

import argparse
import os

from cancer_biomarkers.config import DataConfig
from cancer_biomarkers.pipeline import run


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--data",
        choices=["synthetic", "tcga_brca"],
        default="synthetic",
        help="Which data source to use (default: synthetic).",
    )
    p.add_argument(
        "--nested",
        action="store_true",
        help="Also run nested cross-validation (slower).",
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    if args.nested:
        os.environ["CB_RUN_NESTED_CV"] = "1"
    run(data_cfg=DataConfig(source=args.data))


if __name__ == "__main__":
    main()
