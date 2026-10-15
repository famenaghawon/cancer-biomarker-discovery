"""Feature selection must not depend on hold-out samples."""
import numpy as np
import pandas as pd

from cancer_biomarkers.ml import _select_top_k


def test_select_top_k_ignores_holdout_labels():
    rng = np.random.default_rng(0)
    n, g = 40, 20
    X = pd.DataFrame(rng.normal(size=(n, g)),
                     columns=[f"G{i}" for i in range(g)])
    y = pd.Series([0] * 20 + [1] * 20)

    train_mask = np.zeros(n, dtype=bool)
    train_mask[:30] = True

    # Strong signal only on the "train" half
    X.loc[train_mask & (y == 1), "G0"] += 5.0

    X_train = X.loc[train_mask]
    y_train = y.loc[train_mask]
    top_a = _select_top_k(X_train, y_train, k=3, seed=42)

    # Corrupt hold-out labels: selection is computed on train only, so the
    # result must not change.
    y_corrupt = y.copy()
    y_corrupt.loc[~train_mask] = 1 - y_corrupt.loc[~train_mask]
    top_b = _select_top_k(X_train, y_train, k=3, seed=42)

    assert top_a == top_b
    assert top_a[0] == "G0"


def test_select_top_k_deterministic_with_seed():
    rng = np.random.default_rng(1)
    X = pd.DataFrame(rng.normal(size=(30, 10)),
                     columns=[f"G{i}" for i in range(10)])
    y = pd.Series([0] * 15 + [1] * 15)
    X.loc[y == 1, "G3"] += 3.0

    a = _select_top_k(X, y, k=5, seed=42)
    b = _select_top_k(X, y, k=5, seed=42)
    assert a == b
    assert a[0] == "G3"
