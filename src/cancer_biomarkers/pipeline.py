"""End-to-end orchestration."""
from __future__ import annotations

import json
from dataclasses import asdict

import joblib

from .config import PATHS, DataConfig, DEConfig, MLConfig, PathConfig
from .data import generate_dataset
from .differential_expression import differential_expression
from .logging_utils import get_logger
from .ml import count_true_in_top_k, train_and_evaluate
from .plotting import heatmap_top_genes, roc_curve_plot, volcano_plot

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

    # 3) ML (fast path)
    ml_results = train_and_evaluate(ds.expression, ds.labels, ml_cfg)

    from .ml import nested_cv_evaluate
    nested_results = nested_cv_evaluate(ds.expression, ds.labels, ml_cfg)

    # 3b) Recovery stats for top-k
    recovery = count_true_in_top_k(ml_results["top_genes"], ds.true_biomarkers)

    # 4) Figures
    volcano_plot(
        de_df,
        paths.figures_dir / "volcano.png",
        min_abs_log2fc=de_cfg.min_abs_log2fc,
        alpha=de_cfg.alpha,
    )
    top = ml_results["top_genes"][:25]
    heatmap_top_genes(ds.expression, ds.labels, top, paths.figures_dir / "heatmap.png")

    # ROC curve for the best-performing model on the held-out test set
    # (regenerate a small fit here so we have the test-set probabilities)
    from sklearn.model_selection import train_test_split
    X_tr, X_te, y_tr, y_te = train_test_split(
        ds.expression, ds.labels,
        test_size=ml_cfg.test_size,
        stratify=ds.labels,
        random_state=ml_cfg.random_state,
    )
    from .ml import _make_models, _select_top_k
    genes = _select_top_k(X_tr, y_tr, ml_cfg.top_k_features, ml_cfg.random_state)
    best = _make_models(ml_cfg)["logistic_regression"]
    best.fit(X_tr[genes], y_tr)
    y_proba = best.predict_proba(X_te[genes])[:, 1]
    roc_curve_plot(y_te.values, y_proba,
                   paths.figures_dir / "roc_curve.png",
                   label="logistic_regression")

    # 5) Persist everything
    summary = {
        "config": {
            "data": asdict(data_cfg),
            "de": asdict(de_cfg),
            "ml": asdict(ml_cfg),
        },
        "n_significant_biomarkers": int(de_df["significant"].sum()),
        "true_biomarkers": sorted(ds.true_biomarkers),
        "true_biomarkers_recovered": sorted(
            set(ds.true_biomarkers) & set(de_df.loc[de_df["significant"], "gene"])),
        "top_k_recovery": recovery,
        "ml": ml_results,
        "ml_nested": nested_results,
    }
    with open(paths.results_dir / "analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    joblib.dump(ml_results, paths.results_dir / "ml_results.pkl")

    logger.info(
        "Top-%d recovery: %d/%d true biomarkers | %s",
        recovery["k"], recovery["n_true_in_top_k"],
        len(ds.true_biomarkers), ", ".join(recovery["true_hits"]) or "(none)",
    )
    logger.info("✅ Pipeline complete. Results in %s", paths.results_dir)
    return summary
