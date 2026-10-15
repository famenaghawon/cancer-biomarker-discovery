"""Simple batch-effect correction.

A lightweight, dependency-free ComBat-style adjustment: for each gene,
subtract the batch mean and add the global mean. Not as sophisticated as
``pycombat`` or limma's ``removeBatchEffect``, but sufficient to demonstrate
awareness of the confound and to compare pipelines with/without correction.
"""
from __future__ import annotations

import pandas as pd

from .logging_utils import get_logger

logger = get_logger(__name__)


def correct_batch(
    expression: pd.DataFrame,
    batch: pd.Series,
) -> pd.DataFrame:
    """Remove per-batch mean shifts from each gene.

    Parameters
    ----------
    expression
        samples x genes DataFrame on a log scale.
    batch
        sample-level batch labels, indexed by sample id.

    Returns
    -------
    A new DataFrame with batch means removed.
    """
    batch = batch.reindex(expression.index)
    if batch.isna().any():
        raise ValueError("Batch labels missing for some samples")

    global_mean = expression.mean(axis=0)
    corrected = expression.copy()
    for b in batch.unique():
        mask = batch == b
        batch_mean = expression.loc[mask].mean(axis=0)
        corrected.loc[mask] = expression.loc[mask] - batch_mean + global_mean

    logger.info(
        "Batch correction applied: %d batches, %d samples",
        batch.nunique(), len(expression),
    )
    return corrected
