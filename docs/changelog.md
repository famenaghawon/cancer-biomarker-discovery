# Changelog

## 2026-10-16 — Dashboard and docs site

- Interactive Streamlit dashboard (`dashboard/app.py`).
- GitHub Pages documentation site via MkDocs Material.
- Streamlit 1.40+ `width` API migration.

## 2026-10-16 — Methylation source

- Synthetic methylation cohort with real 450k CpG probe IDs.
- M-value transform so the DE + ML stack runs unchanged on methylation data.
- Silent-tumour and moderate-effect design so classification is not saturated.

## 2026-10-16 — Batch correction

- `--batch-correct` flag on the pipeline CLI.
- ComBat-style per-batch mean-shift removal in `batch.correct_batch`.
- Tests verifying zero between-batch spread after correction.

## 2026-10-16 — Survival analysis

- Cox proportional hazards fit via `lifelines`.
- Kaplan–Meier curves stratified by high/low expression.
- Real BRCA clinical metadata (vital status, days to death/follow-up).

## 2026-10-16 — Multi-source data layer

- `data_source` subpackage with synthetic and TCGA-BRCA backends.
- Real public data via UCSC Xena (no auth, cached locally).
- Backwards-compatible `data.py` shim.

## 2026-10-15 — Reproducibility and CI

- GitHub Actions CI on Python 3.10 / 3.11 / 3.12.
- `tasks.ps1` runner.
- Multi-source `--data` flag.
- Fixed packaging (`pyproject.toml`, `requirements.txt`).

## 2026-10-14 — Initial rebuild

- Replaced flat scripts with a proper `src/cancer_biomarkers/` package.
- Correct log₂FC + Welch t-test + BH-FDR for differential expression.
- Leakage-safe feature selection (train-fold only).
- Stratified cross-validation and nested CV.
- Recovery metrics against planted ground truth.
- Synthetic generator with realistic hardness.
