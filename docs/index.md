# Cancer Biomarker Discovery Pipeline

> An end-to-end computational platform for cancer biomarker discovery on
> **synthetic TCGA-like** and **real TCGA-BRCA** data. Combines rigorous
> differential expression (Welch's t-test + Benjamini–Hochberg FDR) with a
> leakage-safe machine-learning workflow, cross-validated and audited by tests
> and CI.

[![CI](https://github.com/famenaghawon/cancer-biomarker-discovery/actions/workflows/ci.yml/badge.svg)](https://github.com/famenaghawon/cancer-biomarker-discovery/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/famenaghawon/cancer-biomarker-discovery/blob/main/LICENSE)

!!! warning "Scope"
    All data is synthetic or drawn from public TCGA-BRCA; results are for
    methodology demonstration only and are not clinically meaningful.

## Quick start

```bash
git clone https://github.com/famenaghawon/cancer-biomarker-discovery.git
cd cancer-biomarker-discovery

python -m venv .venv
source .venv/bin/activate    # or .venv\Scripts\Activate.ps1 on Windows

pip install -e ".[dev]"

# Synthetic cohort (fast, offline, deterministic)
python scripts/run_pipeline.py

# Real TCGA-BRCA subset (downloads ~40 MB on first run)
python scripts/run_pipeline.py --data tcga_brca

# Methylation multi-omics cohort
python scripts/run_pipeline.py --data methylation
