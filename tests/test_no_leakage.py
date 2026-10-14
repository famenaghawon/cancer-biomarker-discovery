import numpy as np
import pandas as pd

from cancer_biomarkers.ml import _select_top_k


def test_select_top_k_ignores_holdout():
    rng = np.random.default_rng(0)
    X = pd.DataFrame(rng.normal(size=(20, 5)), columns=list("ABCDE"))
    y = pd.Series([0] * 10 + [1] * 10)
    X.loc[y == 1, "A"] += 5
    top = _select_top_k(X, y, k=1, seed=0)
    assert top == ["A"]
