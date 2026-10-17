# Cancer Biomarker Discovery Pipeline

> An end-to-end computational platform for cancer biomarker discovery on
> **synthetic TCGA-like** and **real TCGA-BRCA** data. Combines rigorous
> differential expression (Welch's t-test + Benjamini–Hochberg FDR) with a
> leakage-safe machine-learning workflow, cross-validated Cox survival
> analysis, full ComBat batch correction, and multi-omics support — all audited
> by tests and CI.

[![CI](https://github.com/famenaghawon/cancer-biomarker-discovery/actions/workflows/ci.yml/badge.svg)](https://github.com/famenaghawon/cancer-biomarker-discovery/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23433739.svg)](https://doi.org/10.5281/zenodo.23433739)
[![Docs](https://img.shields.io/badge/docs-live-blue)](https://famenaghawon.github.io/cancer-biomarker-discovery/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/famenaghawon/cancer-biomarker-discovery/blob/main/LICENSE)

!!! warning "Scope"
    All data is synthetic or drawn from public TCGA-BRCA; results are for
    methodology demonstration only and are not clinically meaningful. This
    pipeline is a research tool, not a clinical decision support system.

## Quick start

```bash
git clone https://github.com/famenaghawon/cancer-biomarker-discovery.git
cd cancer-biomarker-discovery

# Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate    # or .venv\Scripts\Activate.ps1 on Windows

# Install the package and dev dependencies
pip install -e ".[dev]"

# Synthetic cohort (fast, offline, deterministic)
python scripts/run_pipeline.py
