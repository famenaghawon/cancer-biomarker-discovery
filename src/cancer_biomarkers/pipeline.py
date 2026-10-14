"""End-to-end orchestration."""
from __future__ import annotations

import json
from dataclasses import asdict

import joblib
import pandas as pd

from .config import DataConfig, DEConfig, MLConfig, PATHS, PathConfig
from .data import generate_dataset
from .differential_expression import differential_expression
from .logging_utils import get_logger
from .ml import train_and_evaluate
from .plotting import heatmap_top_genes, volcano_plot

logger = get_logger(__name__)


def run(
    data_cfg: DataConfig | None = None,
    de_cfg: DEConfig | None = None,
    ml_cfg: MLConfig | None = None,
    paths: PathConfig | None = None,
) -> dict:
    data_cfg = data_cfg or DataConfig()
    de_cfg = de_cfg or DEConfig()
    ml_cfg = ml_cfg or MLConfig()
    paths = paths or PATHS
    paths.ensure()

    # 1) Data
    ds = generate_dataset(data_cfg)
    ds.expression.to_parquet(paths.results_dir / "expression_full.parquet")
    ds.clinical.to_csv(paths.results_dir / "clinical_data.csv", index=False)

    # 2) DE
    de_df = differential_expression(ds.expression, ds.labels, de_cfg)
    de_df.to_csv(paths.results_dir / "de_results.csv", index=False)
    de_df[de_df["significant"]].to_csv(
        paths.results_dir / "significant_biomarkers.csv", index=False)

    # 3) ML
    ml_results = train_and_evaluate(ds.expression, ds.labels, ml_cfg)

    # 4) Figures
    volcano_plot(de_df, paths.figures_dir / "volcano.png")
    top = [g for g in ml_results["top_genes"][:25]]
    heatmap_top_genes(ds.expression, ds.labels, top, paths.figures_dir / "heatmap.png")

    # 5) Persist everything
    summary = {
        "config": {
            "data": asdict(data_cfg),
            "de": asdict(de_cfg),
            "ml": asdict(ml_cfg),
        },
        "n_significant_biomarkers": int(de_df["significant"].sum()),
        "true_biomarkers_recovered": sorted(
            set(ds.true_biomarkers) & set(de_df.loc[de_df["significant"], "gene"])),
        "ml": ml_results,
    }
    with open(paths.results_dir / "analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    joblib.dump(ml_results, paths.results_dir / "ml_results.pkl")

    logger.info("✅ Pipeline complete. Results in %s", paths.results_dir)
    return summary
