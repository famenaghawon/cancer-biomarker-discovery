"""Machine learning with proper cross-validation and no data leakage."""
from __future__ import annotations

import numpy as np
import pandas as pd
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
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .config import MLConfig
from .logging_utils import get_logger

logger = get_logger(__name__)


def _make_models(cfg: MLConfig) -> dict[str, Pipeline]:
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


def _select_top_k(
    expression: pd.DataFrame, labels: pd.Series, k: int, seed: int
) -> list[str]:
    """Rank genes by |Welch t-statistic| computed *only on training data*,
    with a small random jitter so feature selection is not deterministic."""
    from scipy import stats
    rng = np.random.default_rng(seed)
    t_stat, _ = stats.ttest_ind(
        expression.loc[labels == 1].values,
        expression.loc[labels == 0].values,
        axis=0, equal_var=False, nan_policy="omit",
    )
    t_stat = np.nan_to_num(t_stat, nan=0.0)
    jitter = rng.normal(0, 0.05, size=t_stat.shape) * np.abs(t_stat).mean()
    order = np.argsort(-np.abs(t_stat + jitter))
    return list(expression.columns[order[:k]])


def train_and_evaluate(
    expression: pd.DataFrame,
    labels: pd.Series,
    cfg: MLConfig | None = None,
) -> dict:
    cfg = cfg or MLConfig()

    X_train, X_test, y_train, y_test = train_test_split(
        expression, labels,
        test_size=cfg.test_size,
        stratify=labels,
        random_state=cfg.random_state,
    )

    # --- Feature selection on TRAINING data only ---
    top_genes = _select_top_k(X_train, y_train, cfg.top_k_features, cfg.random_state)
    X_train_k = X_train[top_genes]
    X_test_k = X_test[top_genes]
    logger.info("Selected top %d genes on training fold only", len(top_genes))

    results: dict = {"top_genes": top_genes, "models": {}}
    cv = StratifiedKFold(n_splits=cfg.n_splits, shuffle=True,
                         random_state=cfg.random_state)

    for name, model in _make_models(cfg).items():
        model.fit(X_train_k, y_train)
        y_pred = model.predict(X_test_k)
        y_proba = (model.predict_proba(X_test_k)[:, 1]
                   if hasattr(model, "predict_proba") else None)

        cv_scores = cross_val_score(model, X_train_k, y_train, cv=cv, scoring="roc_auc")

        metrics = {
            "accuracy": float(accuracy_score(y_test, y_pred)),
            "precision": float(precision_score(y_test, y_pred, zero_division=0)),
            "recall": float(recall_score(y_test, y_pred, zero_division=0)),
            "f1": float(f1_score(y_test, y_pred, zero_division=0)),
            "roc_auc": float(roc_auc_score(y_test, y_proba)) if y_proba is not None else None,
            "cv_roc_auc_mean": float(cv_scores.mean()),
            "cv_roc_auc_std": float(cv_scores.std()),
            "confusion_matrix": confusion_matrix(y_test, y_pred).tolist(),
        }
        results["models"][name] = metrics
        logger.info("%s | acc=%.3f  f1=%.3f  auc=%.3f  cv_auc=%.3f±%.3f",
                    name, metrics["accuracy"], metrics["f1"],
                    metrics["roc_auc"] or float("nan"),
                    metrics["cv_roc_auc_mean"], metrics["cv_roc_auc_std"])

    return results
