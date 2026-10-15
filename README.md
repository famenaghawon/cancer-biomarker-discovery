# 🧬 Cancer Biomarker Discovery Pipeline

[![CI](https://github.com/famenaghawon/cancer-biomarker-discovery/actions/workflows/ci.yml/badge.svg)](https://github.com/famenaghawon/cancer-biomarker-discovery/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

> An end-to-end computational platform for cancer biomarker discovery on
> **synthetic TCGA-like** data. Combines rigorous differential expression
> (Welch's t-test + Benjamini–Hochberg FDR) with a leakage-safe machine-learning
> workflow, cross-validated and audited by tests and CI.

---

## ⚠️ Disclaimer

This is a **demo and research-education tool**. All data is synthetic, and
results are not clinically meaningful. Nothing in this repository should be
used for diagnosis, treatment selection, or any clinical decision.

---

## ✨ What it does

1. **Generates** a reproducible synthetic expression matrix with known
   ground-truth biomarkers, batch effects, silent tumors, and decoy genes.
2. **Runs differential expression**: per-gene log₂ fold change, Welch's t-test,
   and Benjamini–Hochberg FDR control.
3. **Trains classifiers** with feature selection performed **only on training
   data** — no leakage across the train/test boundary.
4. **Reports** both cross-validated and held-out test metrics, plus a
   top-K recovery rate against the planted biomarker panel.
5. **Emits** a volcano plot, a z-scored heatmap, and a ROC curve, along with a
   JSON summary and CSV tables of results.

---

## 🚀 Quick Start

```bash
git clone https://github.com/famenaghawon/cancer-biomarker-discovery.git
cd cancer-biomarker-discovery

# Create and activate a virtual environment
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
# macOS / Linux:
source .venv/bin/activate

# Install the package and dev dependencies (editable install)
pip install -e ".[dev]"

# Run the full pipeline
python scripts/run_pipeline.py
```

Results land in `results/`:

```
results/
├── analysis_summary.json          # config + DE counts + recovery + ML metrics
├── clinical_data.csv              # sample metadata (label, age, stage, batch)
├── de_results.csv                 # full DE table, sorted by q-value
├── significant_biomarkers.csv     # filtered to significant genes
├── expression_full.parquet        # full samples × genes matrix
├── ml_results.pkl                 # pickled ML result dict
└── figures/
    ├── volcano.png                # log2FC vs -log10(q), thresholds aligned to DEConfig
    ├── heatmap.png                # top-25 biomarkers, z-scored on log2 scale
    └── roc_curve.png              # ROC for logistic regression on the test set
```

---

## ⚡ Task Shortcuts (PowerShell)

Once installed, use the bundled task runner for common actions:

```powershell
.\tasks.ps1 run     # run the full pipeline
.\tasks.ps1 test    # pytest -q
.\tasks.ps1 lint    # ruff check src tests
.\tasks.ps1 clean   # remove __pycache__, test caches, and stale results
.\tasks.ps1 help    # show usage
```

If PowerShell blocks script execution the first time in a session:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

---

## 🧪 Tests

```bash
pytest -q
```

| Test | Validates |
|---|---|
| `test_differential_expression.py` | A known log₂FC signal is detected and ranked first |
| `test_no_leakage.py` | Feature selection does not depend on hold-out labels; deterministic under fixed seed |
| `test_biomarker_recovery.py` | Planted biomarkers are recovered at a reasonable rate under default-like settings |
| `test_pipeline_smoke.py` | The pipeline runs end-to-end and writes all expected outputs |

CI runs the suite and a smoke run of the pipeline on **Python 3.10, 3.11, and 3.12**.
The pipeline's `results/` directory is uploaded as a build artefact on every run.

---

## 📁 Project Structure

```
cancer-biomarker-discovery/
├── .github/workflows/ci.yml       # GitHub Actions: lint + test + smoke run on 3.10–3.12
├── docs/
│   ├── methods.md                 # full methodological write-up
│   ├── figures.md                 # what each plot shows and how to sanity-check
│   └── reviewer_faq.md            # answers to skeptical reviewer questions
├── results/                       # outputs (gitignored; only .gitkeep tracked)
├── scripts/
│   └── run_pipeline.py            # CLI entry point
├── src/cancer_biomarkers/
│   ├── __init__.py
│   ├── config.py                  # dataclasses: DataConfig, DEConfig, MLConfig, PathConfig
│   ├── data.py                    # synthetic TCGA-like data generator
│   ├── differential_expression.py # log2FC + Welch t-test + BH-FDR
│   ├── logging_utils.py           # consistent log formatting
│   ├── ml.py                      # leakage-safe ML, CV, nested CV, recovery
│   ├── pipeline.py                # end-to-end orchestration
│   └── plotting.py                # volcano, heatmap, ROC figures
├── tests/
│   ├── test_biomarker_recovery.py
│   ├── test_differential_expression.py
│   ├── test_no_leakage.py
│   └── test_pipeline_smoke.py
├── tasks.ps1                      # PowerShell task runner
├── pyproject.toml
├── requirements.txt
├── requirements-dev.txt
├── README.md                      # this file
├── PROJECT_SUMMARY.md             # short summary + change history
└── LICENSE
```

---

## 🔬 Methods (short version)

The full write-up is in [`docs/methods.md`](docs/methods.md). A one-paragraph
summary of each stage:

### Synthetic data generation

The generator produces a samples × genes matrix on a **log₂ scale**, then
converts to linear units for storage. Each entry combines:

```
log2_expr[i, g] = baseline[g] + noise[i, g] + lib[i] + batch[i] + tumor_effect[i, g]
```

- `baseline[g] ~ N(8, 1.0)` — typical log₂ expression.
- `noise[i, g] ~ N(0, (noise_sd · per_gene_sd[g])²)` — heterogeneous per-gene variance.
- `lib[i]`, `batch[i]` — library-size shift and 3-batch structure.
- `tumor_effect[i, g]` — applied only to the ~25 true biomarker panel, drawn
  from `N(effect_size, (effect_size/2)²)`, direction flipped 25% of the time.

Realism features: **20% of tumors are "silent"** on the biomarker panel, and
**~25 decoy genes** receive weak label-correlated shifts to force honest
feature selection. All randomness flows through `numpy.random.default_rng(seed)`.

### Differential expression

Per gene *g*:

```
log2FC[g] = log2( (mean_tumor[g] + ε) / (mean_normal[g] + ε) )
```

- **Welch's t-test** (unequal variances) via `scipy.stats.ttest_ind(equal_var=False)`.
- **Benjamini–Hochberg FDR** via `statsmodels.stats.multitest.multipletests(method="fdr_bh")`.
- **Significance rule**: `q_value < alpha` **and** `|log2FC| >= min_abs_log2fc`.
  Defaults: `alpha = 0.05`, `min_abs_log2fc = 0.6`.

### Machine learning

- **Stratified split** — 80/20 train/test, class ratio preserved.
- **Leakage-safe feature selection** — top-K genes by |Welch t-statistic| computed
  **only on the training fold**; a small jitter avoids deterministic selection.
- **Models** — logistic regression (`C=0.1`, StandardScaler) and random forest
  (`n_estimators=200, max_depth=5, min_samples_leaf=5`).
- **Fast path** — one stratified split; 5-fold stratified CV on the training
  features; final metrics on the held-out test set.
- **Strict path** — `nested_cv_evaluate` re-selects features inside each outer
  fold, eliminating residual optimism. Slower; use it for reporting.
- **Recovery metric** — `count_true_in_top_k` reports how many of the selected
  top-K genes match the planted biomarker panel.

### Figures

- **Volcano plot** — log₂FC vs. −log₁₀(q); threshold lines match `DEConfig`.
- **Heatmap** — top-25 biomarkers, log₂-transformed and z-scored per gene.
- **ROC curve** — logistic regression on the held-out test set; the AUC shown
  in the legend matches `analysis_summary.json`.

See [`docs/figures.md`](docs/figures.md) for how to read each figure and how
to sanity-check it against the JSON numbers.

---

## ⚙️ Configuration

All tunables live in `src/cancer_biomarkers/config.py` as frozen dataclasses.
Every value used in a run is serialised into `results/analysis_summary.json`
under `config`, so a report can be reproduced exactly from the summary alone.

| Field | Default | Purpose |
|---|---|---|
| `DataConfig.n_samples` | `200` | Number of samples (tumour + normal) |
| `DataConfig.n_genes` | `500` | Number of genes in the expression matrix |
| `DataConfig.n_true_biomarkers` | `25` | Size of the planted biomarker panel |
| `DataConfig.effect_size` | `0.5` | log₂-scale mean shift for true biomarkers |
| `DataConfig.noise_sd` | `1.0` | Per-gene noise multiplier |
| `DataConfig.batch_effect_sd` | `0.25` | Library-size and batch-effect magnitude |
| `DEConfig.alpha` | `0.05` | FDR significance threshold |
| `DEConfig.min_abs_log2fc` | `0.6` | Minimum |log₂FC| for a significant call |
| `MLConfig.test_size` | `0.2` | Held-out fraction |
| `MLConfig.n_splits` | `5` | Stratified CV folds |
| `MLConfig.top_k_features` | `30` | Number of genes selected for ML |

**Effect-size → expected performance** (documented in `docs/methods.md`):

| `effect_size` | DE hits | ML accuracy | Use |
|---|---|---|---|
| `1.0` | 20–30 | 0.98–1.00 | Sanity check only |
| `0.7` | 10–15 | 0.93–0.96 | Clean, well-powered demo |
| **`0.5`** | **6–10** | **0.83–0.88** | **Default — honest demo** |
| `0.4` | 2–4 | 0.75–0.82 | Adversarially hard |

### Custom runs

```python
from cancer_biomarkers.config import DataConfig, DEConfig, MLConfig
from cancer_biomarkers.pipeline import run

summary = run(
    data_cfg=DataConfig(effect_size=0.7, n_samples=300),
    de_cfg=DEConfig(min_abs_log2fc=0.8),
    ml_cfg=MLConfig(top_k_features=50),
)
```

---

## 🧠 Design principles

1. **No leakage.** Feature selection, scaling, and hyperparameter choices use
   only training-fold data. Enforced by `tests/test_no_leakage.py`.
2. **Honest defaults.** The synthetic data is deliberately not trivially
   separable. A pipeline that only works because the data is easy is worse
   than no pipeline.
3. **Everything configurable, nothing hard-coded.** Thresholds live in
   dataclasses and appear in the output JSON.
4. **Reproducibility.** Same seed ⇒ identical outputs. All randomness uses
   `numpy.random.default_rng`.
5. **Claims backed by tests.** Every methodological claim in the docs has a
   test that would fail if the claim were violated.

---

## 📚 Documentation

- [**Methods**](docs/methods.md) — full methodological write-up with formulas,
  thresholds, and design decisions.
- [**Figures**](docs/figures.md) — what each plot shows, how it is built, how
  to sanity-check it against the summary JSON.
- [**Reviewer FAQ**](docs/reviewer_faq.md) — answers to twelve skeptical
  reviewer questions (why synthetic data? no leakage? why Welch? etc.).
- [**Project Summary**](PROJECT_SUMMARY.md) — high-level overview and change
  history.

---

## 🔬 Research applications

This pipeline is a compact demonstration of the standard workflow used in
biomarker-discovery research. It is directly applicable as a starting point for:

- **Early detection** — evaluating candidate biomarkers for cancer screening.
- **Prognostic stratification** — separating patients by predicted outcome.
- **Therapeutic targeting** — prioritising candidate drug targets from
  expression data.
- **Personalised medicine** — treatment selection based on molecular profiles.

The same code structure (data → DE → ML → figures → summary) works on any
continuous expression matrix. Swapping the synthetic generator for a real
data loader is the primary extension point; see the "Limitations" section of
`docs/methods.md` for what that change would require.

---

## 🛠 Requirements

- **Python**: 3.10, 3.11, or 3.12 (CI-tested). 3.13/3.14 may work but are not
  currently guaranteed by the test matrix.
- **Runtime**: `pandas`, `numpy`, `scipy`, `scikit-learn`, `statsmodels`,
  `matplotlib`, `seaborn`, `joblib`, `pyarrow`.
- **Dev**: `pytest`, `pytest-cov`, `ruff`.

All pinned with lower bounds in `requirements.txt` and `requirements-dev.txt`.
The full dependency set also lives in `pyproject.toml` under `[project]`
and `[project.optional-dependencies].dev`.

---

## 🤝 Contributing

Contributions are welcome. Please:

1. Fork the repository.
2. Create a feature branch: `git checkout -b feature/my-feature`.
3. Ensure `ruff check src tests` and `pytest -q` pass locally.
4. Commit with a descriptive message.
5. Open a pull request.

The CI workflow will run the full test matrix on your PR.

---

## 📝 Citation

If you use this software in research or teaching, please cite:

```bibtex
@software{cancer_biomarker_2025,
  title  = {Cancer Biomarker Discovery Pipeline},
  author = {Famenaghawon, Freedom},
  year   = {2025},
  url    = {https://github.com/famenaghawon/cancer-biomarker-discovery}
}
```

---

## 📧 Contact

**Freedom Famenaghawon** — f.e.amenaghawon@gmail.com

For questions about methods, see [`docs/reviewer_faq.md`](docs/reviewer_faq.md)
first — most reviewer questions are answered there.

---

## 📄 License

MIT — see [`LICENSE`](LICENSE) for the full text.

---

<p align="center">
  <sub>Built as a demonstration of rigorous, reproducible biomarker-discovery methodology.</sub><br>
  <sub>Synthetic data only. Not for clinical use.</sub>
</p>
