"""Machine learning with proper cross-validation and no data leakage."""
from __future__ import annotations

from collections import Counter
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import (
    StratifiedKFold,
    cross_val_score,
    train_test_split,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .config import MLConfig
from .logging_utils import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Model definitions
# ---------------------------------------------------------------------------

def _make_models(cfg: MLConfig) -> dict[str, Pipeline]:
    """Return the baseline classifier pipelines.

    Hyperparameters are intentionally fixed so results are reproducible across
    runs. A hyperparameter search belongs inside the inner loop of nested CV
    (see ``nested_cv_evaluate``); it is not done here so the fast path used by
    the smoke test stays under one second.
    """
    return {
        "logistic_regression": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(
                C=0.1,
                max_iter=2000,
                random_state=cfg.random_state,
            )),
        ]),
        "random_forest": Pipeline([
            ("clf", RandomForestClassifier(
                n_estimators=200,
                max_depth=5,
                min_samples_leaf=5,
                random_state=cfg.random_state,
                n_jobs=1,
            )),
        ]),
    }


# ---------------------------------------------------------------------------
# Feature selection
# ---------------------------------------------------------------------------

def _select_top_k(
    expression: pd.DataFrame,
    labels: pd.Series,
    k: int,
    seed: int,
) -> list[str]:
    """Rank genes by |Welch t-statistic| computed *only on training data*.

    A small random jitter is added to the t-statistics so that feature
    selection has minor variability across runs and cannot be perfectly
    gamed by the classifier. The jitter is scaled by the mean |t| so it
    is meaningful for any dataset size.
    """
    rng = np.random.default_rng(seed)
    t_stat, _ = stats.ttest_ind(
        expression.loc[labels == 1].values,
        expression.loc[labels == 0].values,
        axis=0,
        equal_var=False,
        nan_policy="omit",
    )
    t_stat = np.nan_to_num(t_stat, nan=0.0)
    jitter = rng.normal(0, 0.05, size=t_stat.shape) * np.abs(t_stat).mean()
    order = np.argsort(-np.abs(t_stat + jitter))
    return list(expression.columns[order[:k]])


def count_true_in_top_k(
    top_genes: list[str],
    true_biomarkers: list[str],
) -> dict[str, Any]:
    """How many of the selected top-k genes are known true biomarkers.

    Accuracy tells you whether the classifier works; recovery tells you
    whether it found the right thing. A model can be right for the wrong
    reasons — this metric closes that gap.
    """
    true_set = set(true_biomarkers)
    hits = [g for g in top_genes if g in true_set]
    return {
        "k": len(top_genes),
        "n_true_in_top_k": len(hits),
        "true_hits": hits,
        "recovery_rate": len(hits) / max(1, len(true_set)),
    }


# ---------------------------------------------------------------------------
# Fast path: single split + CV on the training features
# ---------------------------------------------------------------------------

def train_and_evaluate(
    expression: pd.DataFrame,
    labels: pd.Series,
    cfg: MLConfig | None = None,
) -> dict[str, Any]:
    """Fast evaluation: one stratified split, one CV on training features.

    Feature selection is performed on the training fold only. The held-out
    test set is never used for selection, scaling, or model fitting.
    """
    cfg = cfg or MLConfig()

    X_train, X_test, y_train, y_test = train_test_split(
        expression,
        labels,
        test_size=cfg.test_size,
        stratify=labels,
        random_state=cfg.random_state,
    )

    # --- Feature selection on TRAINING data only ---
    top_genes = _select_top_k(
        X_train, y_train, cfg.top_k_features, cfg.random_state
    )
    X_train_k = X_train[top_genes]
    X_test_k = X_test[top_genes]
    logger.info("Selected top %d genes on training fold only", len(top_genes))

    results: dict[str, Any] = {"top_genes": top_genes, "models": {}}
    cv = StratifiedKFold(
        n_splits=cfg.n_splits,
        shuffle=True,
        random_state=cfg.random_state,
    )

    for name, model in _make_models(cfg).items():
        model.fit(X_train_k, y_train)
        y_pred = model.predict(X_test_k)
        y_proba = (
            model.predict_proba(X_test_k)[:, 1]
            if hasattr(model, "predict_proba") else None
        )

        cv_scores = cross_val_score(
            model, X_train_k, y_train, cv=cv, scoring="roc_auc"
        )

        metrics = {
            "accuracy": float(accuracy_score(y_test, y_pred)),
            "precision": float(precision_score(y_test, y_pred, zero_division=0)),
            "recall": float(recall_score(y_test, y_pred, zero_division=0)),
            "f1": float(f1_score(y_test, y_pred, zero_division=0)),
            "roc_auc": float(roc_auc_score(y_test, y_proba))
                        if y_proba is not None else None,
            "cv_roc_auc_mean": float(cv_scores.mean()),
            "cv_roc_auc_std": float(cv_scores.std()),
            "confusion_matrix": confusion_matrix(y_test, y_pred).tolist(),
        }
        results["models"][name] = metrics
        logger.info(
            "%s | acc=%.3f  f1=%.3f  auc=%.3f  cv_auc=%.3f±%.3f",
            name,
            metrics["accuracy"],
            metrics["f1"],
            metrics["roc_auc"] or float("nan"),
            metrics["cv_roc_auc_mean"],
            metrics["cv_roc_auc_std"],
        )

    return results


# ---------------------------------------------------------------------------
# Strict path: nested cross-validation
# ---------------------------------------------------------------------------

def nested_cv_evaluate(
    expression: pd.DataFrame,
    labels: pd.Series,
    cfg: MLConfig | None = None,
) -> dict[str, Any]:
    """Outer CV for performance; feature selection happens inside each fold.

    This avoids the optimistic bias of selecting features once on the full
    training set before cross-validation. It is slower than
    ``train_and_evaluate`` and is intended for deeper reports, not the default
    smoke path.
    """
    cfg = cfg or MLConfig()

    outer = StratifiedKFold(
        n_splits=cfg.n_splits,
        shuffle=True,
        random_state=cfg.random_state,
    )
    models = _make_models(cfg)

    fold_metrics: dict[str, list[dict[str, float]]] = {
        name: [] for name in models
    }
    all_selected: list[list[str]] = []

    for fold_id, (tr, te) in enumerate(outer.split(expression, labels)):
        X_tr, X_te = expression.iloc[tr], expression.iloc[te]
        y_tr, y_te = labels.iloc[tr], labels.iloc[te]

        # Re-select features inside each fold's training split only.
        genes = _select_top_k(
            X_tr, y_tr, cfg.top_k_features, cfg.random_state + fold_id
        )
        all_selected.append(genes)
        X_tr_k, X_te_k = X_tr[genes], X_te[genes]

        for name, model in models.items():
            m = clone(model)  # fresh estimator each fold
            m.fit(X_tr_k, y_tr)
            y_pred = m.predict(X_te_k)
            y_proba = (
                m.predict_proba(X_te_k)[:, 1]
                if hasattr(m, "predict_proba") else None
            )
            fold_metrics[name].append({
                "accuracy": float(accuracy_score(y_te, y_pred)),
                "f1": float(f1_score(y_te, y_pred, zero_division=0)),
                "roc_auc": float(roc_auc_score(y_te, y_proba))
                            if y_proba is not None else float("nan"),
            })

    summary: dict[str, Any] = {"models": {}, "selection_stability": {}}
    for name, folds in fold_metrics.items():
        aucs = [f["roc_auc"] for f in folds]
        summary["models"][name] = {
            "roc_auc_mean": float(np.mean(aucs)),
            "roc_auc_std": float(np.std(aucs)),
            "accuracy_mean": float(np.mean([f["accuracy"] for f in folds])),
            "f1_mean": float(np.mean([f["f1"] for f in folds])),
            "folds": folds,
        }

    # How often each gene was selected across outer folds.
    counts = Counter(g for genes in all_selected for g in genes)
    summary["selection_stability"] = dict(counts.most_common())

    return summary
