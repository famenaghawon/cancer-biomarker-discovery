import numpy as np
import pandas as pd

from cancer_biomarkers.differential_expression import differential_expression


def test_de_detects_known_signal():
    rng = np.random.default_rng(0)
    n, g = 40, 30
    expr = pd.DataFrame(rng.normal(8, 1, size=(n, g)),
                        columns=[f"G{i}" for i in range(g)])
    labels = pd.Series([0] * 20 + [1] * 20)
    # strong signal in G0
    expr.loc[labels == 1, "G0"] += 4.0

    de = differential_expression(expr, labels)
    top = de.iloc[0]
    assert top["gene"] == "G0"
    assert top["log2fc"] > 0
    assert top["q_value"] < 0.05
