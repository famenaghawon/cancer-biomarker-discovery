"""Download and load a curated TCGA-BRCA subset from UCSC Xena.

The loader is network-optional: files are cached under ``data/cache/`` so
subsequent runs are offline. CI does not use this module; local runs do.

Data source
-----------
UCSC Xena Browser, TCGA-BRCA cohort (public, no auth).

* Expression: log2(norm_count+1), one row per gene, one column per sample.
* Clinical:  sample-level metadata including vital_status, days_to_death,
             days_to_last_followup, pathologic_stage, sample_type.

URLs are stable but not guaranteed forever. If a URL breaks, the loader
falls back to alternatives in ``_URL_CANDIDATES`` and raises a clear error
if none respond.
"""
from __future__ import annotations

import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import DataConfig
from ..logging_utils import get_logger

logger = get_logger(__name__)


# Primary URLs; the loader tries them in order.
_EXPR_URLS = [
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz",
]
_CLIN_URLS = [
    "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
]

_CACHE_DIR = Path("data") / "cache"

# A curated cancer-gene panel keeps the pipeline tractable on a laptop.
# Presence in this list is not a claim of clinical relevance.
_CURATED_PANEL = [
    "BRCA1", "BRCA2", "TP53", "EGFR", "ERBB2", "PIK3CA", "PTEN", "RB1",
    "CDH1", "GATA3", "MAP3K1", "AKT1", "MYC", "CCND1", "MDM2", "CDK4",
    "CDK6", "FGFR1", "FGFR2", "FGFR3", "KRAS", "BRAF", "ALK", "MET",
    "RET", "NF1", "SMAD4", "APC", "VHL", "BAX", "BCL2", "MTOR",
    "PIK3R1", "FOXA1", "TBX3", "RUNX1", "CBFB", "ARID1A", "ARID1B",
    "KMT2C", "NCOR1", "SPEN", "CTCFF", "ZNF217", "AURKA", "AURKB",
    "BUB1", "BUB1B", "MAD2L1", "PLK1", "TTK", "NDC80", "ESR1", "PGR",
    "PRLR", "CYP19A1", "HSD17B1", "IGF1R", "INSR", "IRS1", "IRS2",
]


@dataclass
class TCGADataset:
    expression: pd.DataFrame           # samples x genes, log2 scale
    labels: pd.Series                  # 0 = normal-adjacent, 1 = primary tumor
    clinical: pd.DataFrame             # sample-level metadata
    true_biomarkers: list[str]         # curated panel intersected with genes
    dataset_name: str = "TCGA-BRCA"


def _download(url: str, dest: Path, timeout: int = 60) -> Path:
    """Download ``url`` to ``dest`` unless ``dest`` already exists."""
    if dest.exists() and dest.stat().st_size > 0:
        logger.info("Using cached %s", dest)
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading %s -> %s", url, dest)
    req = urllib.request.Request(url, headers={"User-Agent": "cancer-biomarkers/0.3"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = resp.read()
    if not data:
        raise RuntimeError(f"Empty response from {url}")
    dest.write_bytes(data)
    logger.info("Saved %d bytes", len(data))
    return dest


def _try_download_first(urls: list[str], dest: Path) -> Path:
    last_err: Exception | None = None
    for url in urls:
        try:
            return _download(url, dest)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Download failed for %s: %s", url, exc)
            last_err = exc
    raise RuntimeError(f"All download URLs failed for {dest.name}") from last_err


def _load_expression(cache: Path) -> pd.DataFrame:
    """Return a samples x genes DataFrame, log2 scale, symbols as columns."""
    path = cache / "tcga_brca_expression.tsv.gz"
    _try_download_first(_EXPR_URLS, path)

    df = pd.read_csv(path, sep="\t", index_col=0, compression="infer")
    # Xena files are genes x samples; transpose to samples x genes
    df = df.T
    df.index.name = "sample_id"
    # Drop the first row if it is a duplicate header (Xena quirk)
    if df.index[0] == "sampleID":
        df = df.iloc[1:]
    return df


def _load_clinical(cache: Path) -> pd.DataFrame:
    path = cache / "tcga_brca_clinical.tsv"
    _try_download_first(_CLIN_URLS, path)
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index.name = "sample_id"
    return df


def _build_labels(clinical: pd.DataFrame) -> pd.Series:
    """1 = primary tumor, 0 = solid tissue normal. Others dropped."""
    st = clinical["sample_type"].str.lower()
    is_tumor = st.str.contains("primary tumor", na=False)
    is_normal = st.str.contains("solid tissue normal", na=False)
    keep = is_tumor | is_normal
    labels = pd.Series(
        np.where(is_tumor[keep], 1, 0),
        index=clinical.index[keep],
        name="label",
    )
    return labels


def _build_clinical_summary(clinical: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    """Extract the subset of clinical columns relevant to our pipeline."""
    keep_cols = [
        "sample_type", "pathologic_stage", "vital_status",
        "days_to_death", "days_to_last_followup",
        "age_at_initial_pathologic_diagnosis",
        "gender",
    ]
    present = [c for c in keep_cols if c in clinical.columns]
    sub = clinical.loc[labels.index, present].copy()

    # Derive a clean event/time pair for the survival module.
    # TCGA / Xena export uses "DECEASED" / "LIVING"; other mirrors may use
    # "Dead" / "Alive" or "deceased" / "living". Normalize and accept a
    # superset of spellings so the loader is robust across sources.
    vs_raw = sub.get("vital_status", pd.Series(index=sub.index, dtype=object))
    vs = vs_raw.astype(str).str.lower().str.strip()
    deceased = vs.isin({"dead", "deceased"})
    sub["event"] = deceased.astype(int)

    dtd = pd.to_numeric(sub.get("days_to_death"), errors="coerce")
    dtf = pd.to_numeric(sub.get("days_to_last_followup"), errors="coerce")
    # Event samples use days_to_death; censored samples use last follow-up.
    sub["time_days"] = dtd.where(sub["event"] == 1, dtf)

    # Simplify stage to a Roman numeral
    if "pathologic_stage" in sub.columns:
        sub["stage"] = (
            sub["pathologic_stage"]
            .astype(str)
            .str.extract(r"Stage\s+([IVX]+)", expand=False)
            .fillna("NA")
        )
    else:
        sub["stage"] = "NA"

    sub["label"] = labels
    return sub


def load_tcga_brca(
    cfg: DataConfig | None = None,
    gene_panel: list[str] | None = None,
) -> TCGADataset:
    """Load TCGA-BRCA expression + clinical, filtered to a gene panel.

    Parameters
    ----------
    cfg
        Unused for now (kept for API symmetry with ``generate_dataset``).
    gene_panel
        Genes to keep. Defaults to :data:`_CURATED_PANEL`.

    Returns
    -------
    TCGADataset
    """
    _ = cfg  # reserved for future filtering options
    gene_panel = gene_panel or _CURATED_PANEL

    cache = _CACHE_DIR
    cache.mkdir(parents=True, exist_ok=True)

    expr = _load_expression(cache)
    clinical = _load_clinical(cache)

    labels = _build_labels(clinical)
    # Intersect samples present in both
    common = expr.index.intersection(labels.index)
    expr = expr.loc[common]
    labels = labels.loc[common]

    # Filter to genes actually present
    genes = [g for g in gene_panel if g in expr.columns]
    if not genes:
        raise RuntimeError(
            "None of the curated gene panel is present in the expression matrix"
        )
    expr = expr[genes].astype(float)

    clinical_summary = _build_clinical_summary(clinical, labels)

    logger.info(
        "Loaded TCGA-BRCA: %d samples (%d tumor / %d normal), %d genes",
        len(labels), int(labels.sum()), int((labels == 0).sum()), len(genes),
    )

    return TCGADataset(
        expression=expr,
        labels=labels,
        clinical=clinical_summary,
        true_biomarkers=genes,
    )
