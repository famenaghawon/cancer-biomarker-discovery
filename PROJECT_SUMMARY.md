# Project Summary

**Goal:** demonstrate a rigorous, reproducible biomarker-discovery workflow on synthetic
TCGA-like data — from data generation to differential expression to ML classification —
without the typical methodological shortcuts (mean-difference "fold change", missing
multiple-testing correction, feature-selection leakage).

**Deliverables**
- `src/cancer_biomarkers/` — installable Python package
- `scripts/run_pipeline.py` — one-command reproduction
- `results/` — CSV/JSON figures from the last run
- `tests/` — pytest suite covering DE, ML leakage, and end-to-end smoke

**What was fixed vs. the original repo**
1. `requirements.txt` is now a valid package list.
2. `ml_pipeline*.py` merged into a single, complete pipeline that writes all outputs.
3. DE now uses log₂FC + Welch t-test + BH-FDR (was difference-of-means).
4. Feature selection now happens **inside** the training fold (was full-data leakage).
5. README no longer advertises files/classes that do not exist.
6. `PROJECT_SUMMARY.md` no longer contains PowerShell heredocs.
7. Added `.gitignore`, `pyproject.toml`, tests, structured logging, type hints.
8. Synthetic data no longer trivially separable (smaller effect + noise + batch effect).

## Later additions (second round)

- Nested cross-validation path (`ml.nested_cv_evaluate`, opt-in via `CB_RUN_NESTED_CV`)
- Top-K biomarker recovery metric (`ml.count_true_in_top_k`) — reported in the
  summary JSON and logged at the end of the pipeline
- ROC curve figure (`figures/roc_curve.png`), driven by predictions returned
  from `train_and_evaluate` — no re-fitting in `pipeline.run`
- Heatmap now log₂-transformed and z-scored per gene
- Volcano plot thresholds read from `DEConfig` instead of hard-coded
- Config-driven noise (`DataConfig.noise_sd`, `DataConfig.batch_effect_sd`)
  now actually wired into `data.generate_dataset`
- Clinical stage is `"NA"` for normal samples
- `n_true_biomarkers` clamped to the pool of real gene symbols with a warning
- Professional contact email unified across README, pyproject, LICENSE, docs
