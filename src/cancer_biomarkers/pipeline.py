"""End-to-end orchestration."""
from __future__ import annotations

import json
import os
from dataclasses import asdict

import joblib

from .config import PATHS, DataConfig, DEConfig, MLConfig, PathConfig
from .data_source import generate_dataset, load_tcga_brca
from .differential_expression import differential_expression
from .logging_utils import get_logger
from .ml import count_true_in_top_k, nested_cv_evaluate, train_and_evaluate
from .plotting import (
    heatmap_top_genes,
    kaplan_meier_plot,
    roc_curve_plot,
    volcano_plot,
)

logger = get_logger(__name__)


# Genes commonly reported as prognostic in BRCA. PGR is intentionally
# excluded: it is missing on a large fraction of TCGA-BRCA samples, which
# destabilises the Cox Hessian even after dropna.
_SURVIVAL_PANEL_CANDIDATES = ["ESR1", "ERBB2", "MKI67"]


def _nested_enabled() -> bool:
    """Nested CV is opt-in to keep the default smoke path fast.

    Enable by setting the ``CB_RUN_NESTED_CV`` environment variable to
    ``1``, ``true``, or ``yes`` before running the pipeline.
    """
    return os.environ.get("CB_RUN_NESTED_CV", "").strip().lower() in {
        "1", "true", "yes", "y",
    }


def _load_dataset(cfg: DataConfig):
    """Dispatch to the configured data source.

    Returns an object exposing ``expression``, ``labels``, ``clinical``,
    and a ``has_planted_biomarkers`` flag. Sources with a planted panel
    additionally expose ``true_biomarkers``; real TCGA-BRCA instead exposes
    a ``gene_panel`` (which is not a planted signal). Methylation data is
    transformed to M-values so the downstream DE + ML pipeline sees the
    same kind of matrix as expression.
    """
    if cfg.source == "synthetic":
        return generate_dataset(cfg)
    if cfg.source == "tcga_brca":
        panel = list(cfg.gene_panel) if cfg.gene_panel else None
        return load_tcga_brca(cfg, gene_panel=panel)
    if cfg.source == "methylation":
        from .data_source.synthetic import SyntheticDataset
        from .multiomics import generate_methylation, to_m_values

        methyl = generate_methylation(cfg)
        m_values = to_m_values(methyl.beta)
        logger.info(
            "Applied M-value transform: %d x %d matrix on log2-odds scale",
            m_values.shape[0], m_values.shape[1],
        )
        return SyntheticDataset(
            expression=m_values,
            labels=methyl.labels,
            clinical=methyl.clinical,
            true_biomarkers=methyl.true_biomarkers,
            has_planted_biomarkers=True,
        )
    raise ValueError(f"Unknown data source: {cfg.source!r}")


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
    ds = _load_dataset(data_cfg)
    expression = ds.expression

    # 1b) Optional: batch correction
    if data_cfg.batch_correct and "batch" in ds.clinical.columns:
        from .batch import correct_batch

        logger.info("Applying batch correction (per-batch mean shift removal)")
        expression = correct_batch(expression, ds.clinical["batch"])
        expression.to_parquet(
            paths.results_dir / "expression_full_batch_corrected.parquet",
        )
    elif data_cfg.batch_correct:
        logger.warning(
            "batch_correct=True but no 'batch' column in clinical data; "
            "skipping correction",
        )

    expression.to_parquet(paths.results_dir / "expression_full.parquet")
    ds.clinical.to_csv(paths.results_dir / "clinical_data.csv", index=False)

    # 2) DE — use the (possibly corrected) local `expression`
    de_df = differential_expression(expression, ds.labels, de_cfg)
    de_df.to_csv(paths.results_dir / "de_results.csv", index=False)
    de_df[de_df["significant"]].to_csv(
        paths.results_dir / "significant_biomarkers.csv", index=False)

    # 3) ML (fast path)
    ml_results = train_and_evaluate(
        expression, ds.labels, ml_cfg, return_test_predictions=True,
    )

    # 3b) Recovery — only meaningful for sources with a planted panel
    if getattr(ds, "has_planted_biomarkers", False):
        recovery = count_true_in_top_k(
            ml_results["top_genes"], ds.true_biomarkers,
        )
    else:
        recovery = None
        logger.info(
            "Recovery metric not applicable: %s has no planted biomarker panel",
            data_cfg.source,
        )

    # 3c) Optional: nested CV (opt-in via CB_RUN_NESTED_CV=1)
    nested_results = None
    if _nested_enabled():
        logger.info("CB_RUN_NESTED_CV set — running nested cross-validation")
        nested_results = nested_cv_evaluate(expression, ds.labels, ml_cfg)

    # 3d) Survival analysis — only when the cohort has time/event columns
    survival_summary = None
    if "time_days" in ds.clinical.columns and "event" in ds.clinical.columns:
        from .survival import cox_ph, kaplan_meier, stratify_by_expression

        preferred = [
            g for g in _SURVIVAL_PANEL_CANDIDATES
            if g in expression.columns
        ]
        panel = preferred if len(preferred) >= 2 else ml_results["top_genes"][:3]
        top_gene = panel[0]

        try:
            cox = cox_ph(
                expression,
                ds.clinical["time_days"],
                ds.clinical["event"],
                top_genes=panel,
            )
            cox.coefficients.to_csv(paths.results_dir / "cox_coefficients.csv")
            stratum = stratify_by_expression(expression[top_gene])
            curves = kaplan_meier(
                ds.clinical["time_days"],
                ds.clinical["event"],
                stratify_by=stratum,
            )
            kaplan_meier_plot(
                curves,
                paths.figures_dir / "kaplan_meier.png",
                title=f"KM by {top_gene} expression",
            )
            survival_summary = {
                "concordance": cox.concordance,
                "panel": panel,
                "top_gene": top_gene,
            }
        except Exception as exc:
            logger.warning("Survival analysis skipped: %s", exc, exc_info=True)

    # 4) Figures
    volcano_plot(
        de_df,
        paths.figures_dir / "volcano.png",
        min_abs_log2fc=de_cfg.min_abs_log2fc,
        alpha=de_cfg.alpha,
    )

    top = ml_results["top_genes"][:25]
    heatmap_top_genes(
        expression, ds.labels, top, paths.figures_dir / "heatmap.png",
    )

    best_name = max(
        ml_results["models"],
        key=lambda n: ml_results["models"][n]["roc_auc"] or 0.0,
    )
    preds = ml_results["test_predictions"][best_name]
    roc_curve_plot(
        preds["y_true"], preds["y_proba"],
        paths.figures_dir / "roc_curve.png",
        label=best_name,
    )

    # 5) Summary — gate planted-panel fields on the flag
    if getattr(ds, "has_planted_biomarkers", False):
        true_biomarkers = sorted(ds.true_biomarkers)
        true_biomarkers_recovered = sorted(
            set(ds.true_biomarkers)
            & set(de_df.loc[de_df["significant"], "gene"])
        )
    else:
        true_biomarkers = []
        true_biomarkers_recovered = []

    summary = {
        "config": {
            "data": asdict(data_cfg),
            "de": asdict(de_cfg),
            "ml": asdict(ml_cfg),
        },
        "n_significant_biomarkers": int(de_df["significant"].sum()),
        "true_biomarkers": true_biomarkers,
        "true_biomarkers_recovered": true_biomarkers_recovered,
        "ml": {k: v for k, v in ml_results.items() if k != "test_predictions"},
    }
    if recovery is not None:
        summary["top_k_recovery"] = recovery
    if nested_results is not None:
        summary["ml_nested"] = nested_results
    if survival_summary is not None:
        summary["survival"] = survival_summary

    with open(paths.results_dir / "analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    joblib.dump(ml_results, paths.results_dir / "ml_results.pkl")

    if recovery is not None:
        logger.info(
            "Top-%d recovery: %d/%d true biomarkers | %s",
            recovery["k"], recovery["n_true_in_top_k"],
            len(ds.true_biomarkers), ", ".join(recovery["true_hits"]) or "(none)",
        )
    logger.info("Best model for ROC: %s", best_name)
    logger.info("✅ Pipeline complete. Results in %s", paths.results_dir)
    return summary
