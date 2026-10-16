"""Batch correction must equalize per-batch means."""
import numpy as np
import pandas as pd
import pytest

from cancer_biomarkers.batch import correct_batch


def test_correct_batch_equalizes_batch_means():
    rng = np.random.default_rng(0)
    n_per_batch = 30
    n_genes = 10
    batches = ["A", "B", "C"]

    rows = []
    batch_labels = []
    for i, b in enumerate(batches):
        offset = np.full(n_genes, i * 2.0)  # batch B +2, batch C +4
        rows.append(rng.normal(0, 1, size=(n_per_batch, n_genes)) + offset)
        batch_labels.extend([b] * n_per_batch)

    expr = pd.DataFrame(
        np.vstack(rows), columns=[f"G{i}" for i in range(n_genes)],
    )
    batch = pd.Series(batch_labels, index=expr.index)

    # Sanity: batches start with very different means
    before = expr.groupby(batch).mean()
    assert (before.max() - before.min()).max() > 2.0

    corrected = correct_batch(expr, batch)

    # After correction, every batch should have the SAME per-gene mean.
    # "Same" means "the between-batch spread is now zero," not "the mean is zero."
    after = corrected.groupby(batch).mean()
    spread = (after.max() - after.min()).abs().max()
    assert spread < 1e-8, f"Between-batch spread after correction: {spread}"


def test_correct_batch_preserves_global_mean():
    rng = np.random.default_rng(1)
    n_per_batch = 20
    n_genes = 5
    rows, labels = [], []
    for i, b in enumerate(["A", "B", "C"]):
        rows.append(rng.normal(0, 1, size=(n_per_batch, n_genes)) + i)
        labels.extend([b] * n_per_batch)

    expr = pd.DataFrame(np.vstack(rows), columns=[f"G{i}" for i in range(n_genes)])
    batch = pd.Series(labels, index=expr.index)

    global_before = expr.mean(axis=0)
    corrected = correct_batch(expr, batch)
    global_after = corrected.mean(axis=0)

    # Batch correction shifts each batch, but the overall mean should be unchanged
    assert np.allclose(global_before, global_after, atol=1e-8)


def test_correct_batch_rejects_missing_labels():
    expr = pd.DataFrame(np.zeros((10, 3)), columns=["A", "B", "C"])
    batch = pd.Series(["X"] * 5, index=expr.index[:5])
    with pytest.raises(ValueError):
        correct_batch(expr, batch)
