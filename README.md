# 🧬 Cancer Biomarker Discovery Pipeline

[![CI](https://github.com/famenaghawon/cancer-biomarker-discovery/actions/workflows/ci.yml/badge.svg)](https://github.com/famenaghawon/cancer-biomarker-discovery/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-live-blue)](https://famenaghawon.github.io/cancer-biomarker-discovery/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

> An end-to-end computational platform for cancer biomarker discovery on
> **synthetic TCGA-like** and **real TCGA-BRCA** data. Combines rigorous
> differential expression (Welch's t-test + Benjamini–Hochberg FDR) with a
> leakage-safe machine-learning workflow, cross-validated Cox survival
> analysis, optional batch correction, and multi-omics support — all audited
> by tests and CI.

---

## ⚠️ Disclaimer

This is a **research-education tool**. Data is either synthetic or drawn from
public TCGA-BRCA resources via the UCSC Xena Browser. Results are for
methodology demonstration and are **not clinically meaningful**. Nothing in
this repository should be used for diagnosis, treatment selection, or any
clinical decision.

---

## ✨ What it does

1. **Loads** one of three data sources: synthetic expression, real
   **TCGA-BRCA**, or synthetic **methylation** with real Illumina 450k CpG
   probe IDs.
2. **Runs differential expression** — per-gene log₂ fold change, Welch's
   t-test, and Benjamini–Hochberg FDR control.
3. **Trains classifiers** with feature selection performed **only on training
   data** — no leakage across the train/test boundary.
4. **Reports** both cross-validated and held-out test metrics, plus a
   top-K recovery rate against the planted biomarker panel (recovery is
   reported for `synthetic` and `methylation` sources, which have a known
   ground truth; it is omitted for real `tcga_brca` data).
5. **Fits a Cox proportional-hazards model** (via `lifelines`) and produces
   Kaplan–Meier curves when the cohort has survival data.
6. **Optionally corrects batch effects** (ComBat-style location-shift removal).
7. **Emits** a volcano plot, a z-scored heatmap, a ROC curve, a Kaplan–Meier
   plot, a JSON summary, and CSV tables of results.
8. **Serves** an interactive Streamlit dashboard and a searchable MkDocs site
   from the same results.

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

# Run the pipeline on the synthetic cohort (fast, offline, deterministic)
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
├── cox_coefficients.csv           # (TCGA only) Cox PH coefficients
└── figures/
    ├── volcano.png                # log2FC vs -log10(q), thresholds aligned to DEConfig
    ├── heatmap.png                # top-25 biomarkers, z-scored on log2 scale
    ├── roc_curve.png              # ROC for the best classifier on the test set
    └── kaplan_meier.png           # (TCGA only) KM by top-gene expression
```

### Run against real TCGA-BRCA data

```bash
# Downloads ~40 MB of public expression + clinical data on the first run,
# then caches under data/cache/. No authentication required (UCSC Xena).
python scripts/run_pipeline.py --data tcga_brca
```

Expected on the real cohort:

```
Loaded TCGA-BRCA: 1211 samples (1097 tumor / 114 normal), 59 genes
DE done: 7/59 significant (q<0.050, |log2FC|>=0.60)
logistic_regression | acc=0.971  f1=0.984  auc=0.995
random_forest       | acc=0.963  f1=0.980  auc=0.988
Cox PH: 1193 samples, 195 events, 2 covariates, concordance=0.525
Top-30 recovery: 30/59 true biomarkers
```

### Run on methylation data (multi-omics)

```bash
python scripts/run_pipeline.py --data methylation
```

Generates a synthetic β-value matrix with real 450k CpG IDs, applies the
M-value transform (`M = log2((β + ε) / (1 − β + ε))`), then runs the same
DE + ML stack.

### Apply batch correction

```bash
python scripts/run_pipeline.py --batch-correct
```

Removes per-batch mean shifts before DE. Only fires when the loaded cohort
carries a `batch` column.

### Enable nested cross-validation

```bash
# Bash / macOS / Linux:
CB_RUN_NESTED_CV=1 python scripts/run_pipeline.py
# Windows PowerShell:
$env:CB_RUN_NESTED_CV = "1"; python scripts/run_pipeline.py
```

Re-selects features inside each outer CV fold, eliminating residual optimism
in the reported performance. Slower; use it for publication-style claims.

---

## 📊 Interactive Dashboard

Once a pipeline run has produced results in `results/`, launch the dashboard:

```bash
pip install -e ".[dashboard]"
streamlit run dashboard/app.py```

The dashboard opens at `http://localhost:8501` and shows:

- Header metrics: data source, significant biomarkers, top-K recovery, best ROC-AUC
- Cox PH concordance and survival panel (when survival data is present)
- Volcano, heatmap, ROC, and Kaplan–Meier figures
- Interactive volcano with hover info (Plotly)
- Sortable DE table with CSV download
- Model comparison bar chart
- Cox PH coefficients table

It is read-only over `results/` and does not re-run any pipeline step.

---

## 📖 Hosted Documentation

The full docs are published at
**[famenaghawon.github.io/cancer-biomarker-discovery](https://famenaghawon.github.io/cancer-biomarker-discovery/)**,
built with MkDocs Material and deployed automatically by GitHub Actions on
every push to `docs/` or `mkdocs.yml`.

Locally:

```bash
pip install -r requirements-docs.txt
mkdocs serve     # serves at http://127.0.0.1:8000
```

The site includes:

- **[Home](https://famenaghawon.github.io/cancer-biomarker-discovery/)** — landing page with badges, Quick Start, and links
- **[Methods](https://famenaghawon.github.io/cancer-biomarker-discovery/methods/)** — full methodological write-up
- **[Figures](https://famenaghawon.github.io/cancer-biomarker-discovery/figures/)** — what each plot shows
- **[Reviewer FAQ](https://famenaghawon.github.io/cancer-biomarker-discovery/reviewer_faq/)** — answers to skeptical reviewer questions
- **[API Reference](https://famenaghawon.github.io/cancer-biomarker-discovery/api/)** — public classes and functions
- **[Changelog](https://famenaghawon.github.io/cancer-biomarker-discovery/changelog/)** — release history

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
| `test_batch.py` | Batch correction zeroes between-batch spread and preserves the global mean |
| `test_methylation.py` | β-values stay in [0, 1]; probe IDs are real 450k IDs; M-value transform inverts correctly |
| `test_pipeline_smoke.py` | The pipeline runs end-to-end, writes all expected outputs, and adds `ml_nested` when `CB_RUN_NESTED_CV=1` |

CI runs the suite and a smoke run of the pipeline on **Python 3.10, 3.11, and
3.12**. The pipeline's `results/` directory is uploaded as a build artefact on
every run.

---

## 📁 Project Structure

```
cancer-biomarker-discovery/
├── .github/workflows/
│   ├── ci.yml                     # lint + test + smoke run on 3.10–3.12
│   └── docs.yml                   # build MkDocs site and deploy to gh-pages
├── .streamlit/config.toml         # Streamlit theme and telemetry opt-out
├── dashboard/
│   └── app.py                     # interactive Streamlit dashboard
├── docs/
│   ├── index.md                   # MkDocs landing page
│   ├── methods.md                 # full methodological write-up
│   ├── figures.md                 # what each plot shows and how to sanity-check
│   ├── reviewer_faq.md            # answers to skeptical reviewer questions
│   ├── api.md                     # public API reference
│   └── changelog.md               # release history
├── results/                       # outputs (gitignored; only .gitkeep tracked)
├── scripts/
│   └── run_pipeline.py            # CLI entry point
├── src/cancer_biomarkers/
│   ├── __init__.py
│   ├── config.py                  # dataclasses: DataConfig, DEConfig, MLConfig, PathConfig
│   ├── data.py                    # backwards-compatible shim for data_source
│   ├── data_source/
│   │   ├── __init__.py
│   │   ├── synthetic.py           # synthetic TCGA-like expression generator
│   │   ├── tcga.py                # real TCGA-BRCA loader (UCSC Xena)
│   │   └── cpg_panel.py           # curated 450k CpG probe IDs
│   ├── differential_expression.py # log2FC + Welch t-test + BH-FDR
│   ├── logging_utils.py           # consistent log formatting
│   ├── ml.py                      # leakage-safe ML, CV, nested CV, recovery
│   ├── survival.py                # Cox PH (lifelines) + Kaplan–Meier
│   ├── batch.py                   # ComBat-style per-batch mean-shift correction
│   ├── multiomics.py              # methylation generator + M-value transform
│   ├── pipeline.py                # end-to-end orchestration
│   └── plotting.py                # volcano, heatmap, ROC, Kaplan–Meier figures
├── tests/
│   ├── test_batch.py
│   ├── test_biomarker_recovery.py
│   ├── test_differential_expression.py
│   ├── test_methylation.py
│   ├── test_no_leakage.py
│   └── test_pipeline_smoke.py
├── tasks.ps1                      # PowerShell task runner
├── mkdocs.yml                     # MkDocs Material site config
├── pyproject.toml
├── requirements.txt
├── requirements-dev.txt
├── requirements-docs.txt
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

### Survival analysis

- **Cox proportional hazards** via `lifelines` on standardized covariates.
  Automatically drops collinear and zero-variance covariates before fitting.
- **Kaplan–Meier curves** computed by hand, stratified by high/low expression
  of the top gene.
- Only runs when the loaded cohort has `time_days` and `event` columns.
  The TCGA-BRCA loader populates these from `days_to_death`,
  `days_to_last_followup`, and `vital_status`.

### Batch-effect correction

When the cohort carries a `batch` column, an optional ComBat-style location
adjustment is applied before DE:

```
log2_expr_corrected[i, g] = log2_expr[i, g]
                            - mean(log2_expr[batch(i), g])
                            + global_mean[g]
```

Enable with `--batch-correct`. The corrected matrix is written alongside the
raw one as `expression_full_batch_corrected.parquet`.

### Multi-omics (methylation)

`--data methylation` loads a synthetic β-value matrix with real Illumina 450k
probe IDs. Values are converted to M-values via:

```
M = log2( (β + ε) / (1 − β + ε) )
```

with ε = 1e-6. The transform puts β ∈ [0, 1] onto an unbounded log-odds scale
that the DE module can process identically to expression data.

### Figures

- **Volcano plot** — log₂FC vs. −log₁₀(q); threshold lines match `DEConfig`.
- **Heatmap** — top-25 biomarkers, log₂-transformed and z-scored per gene.
- **ROC curve** — best classifier on the held-out test set.
- **Kaplan–Meier** — survival by high/low expression of the top gene.

See [`docs/figures.md`](docs/figures.md) for how to read each figure and how
to sanity-check it against the JSON numbers.

---

## ⚙️ Configuration

All tunables live in `src/cancer_biomarkers/config.py` as frozen dataclasses.
Every value used in a run is serialised into `results/analysis_summary.json`
under `config`, so a report can be reproduced exactly from the summary alone.

| Field | Default | Purpose |
|---|---|---|
| `DataConfig.source` | `"synthetic"` | One of `"synthetic"`, `"tcga_brca"`, `"methylation"` |
| `DataConfig.batch_correct` | `False` | Apply per-batch mean-shift removal before DE |
| `DataConfig.n_samples` | `200` | Number of samples (tumour + normal) |
| `DataConfig.n_genes` | `500` | Number of genes (or CpG sites) |
| `DataConfig.tumor_ratio` | `0.5` | Fraction of tumour samples |
| `DataConfig.n_true_biomarkers` | `25` | Size of the planted biomarker panel |
| `DataConfig.effect_size` | `0.5` | log₂-scale mean shift for true biomarkers |
| `DataConfig.noise_sd` | `1.0` | Per-gene noise multiplier |
| `DataConfig.batch_effect_sd` | `0.25` | Library-size and batch-effect magnitude |
| `DataConfig.random_state` | `42` | Seed for reproducibility |
| `DEConfig.alpha` | `0.05` | FDR significance threshold |
| `DEConfig.min_abs_log2fc` | `0.6` | Minimum \|log₂FC\| for a significant call |
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

**Recovery metric scope.** The `top_k_recovery` field in
`analysis_summary.json` reports how many of the selected top-K features match
the planted biomarker panel. It is only meaningful for the `synthetic` and
`methylation` sources — those generators plant a known panel. For real
`tcga_brca` data, the field is omitted because the curated gene list that gets
loaded is not a signal. See
[`docs/methods.md`](docs/methods.md) for the full explanation.

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
2. **Honest defaults.** Synthetic data is deliberately not trivially
   separable. A pipeline that only works because the data is easy is worse
   than no pipeline.
3. **Everything configurable, nothing hard-coded.** Thresholds live in
   dataclasses and appear in the output JSON.
4. **Reproducibility.** Same seed ⇒ identical outputs. All randomness uses
   `numpy.random.default_rng`.
5. **Claims backed by tests.** Every methodological claim in the docs has a
   test that would fail if the claim were violated.
6. **Real data is a first-class path.** TCGA-BRCA is not an afterthought — it
   runs through the same DE, ML, and survival stack as the synthetic cohort.

---

## 📚 Documentation

- [**Documentation Site**](https://famenaghawon.github.io/cancer-biomarker-discovery/) — the same docs, hosted and searchable
- [**Methods**](docs/methods.md) — full methodological write-up with formulas, thresholds, and design decisions
- [**Figures**](docs/figures.md) — what each plot shows, how it is built, how to sanity-check it against the summary JSON
- [**Reviewer FAQ**](docs/reviewer_faq.md) — answers to skeptical reviewer questions
- [**API Reference**](docs/api.md) — public classes and functions
- [**Changelog**](docs/changelog.md) — release history
- [**Project Summary**](PROJECT_SUMMARY.md) — high-level overview

---

## 🔬 Research applications

This pipeline is a compact demonstration of the standard workflow used in
biomarker-discovery research. It is directly applicable as a starting point for:

- **Early detection** — evaluating candidate biomarkers for cancer screening.
- **Prognostic stratification** — separating patients by predicted survival.
- **Therapeutic targeting** — prioritising candidate drug targets from
  expression or methylation data.
- **Personalised medicine** — treatment selection based on molecular profiles.

The same code structure (data → DE → ML → survival → figures → summary) works
on any continuous omics matrix. Adding a new modality requires implementing a
loader that returns a `samples × features` DataFrame plus labels — everything
downstream is reused.

---

## 🛠 Requirements

- **Python**: 3.10, 3.11, or 3.12 (CI-tested). 3.13/3.14 may work but are not
  currently guaranteed by the test matrix.
- **Runtime**: `pandas`, `numpy`, `scipy`, `scikit-learn`, `statsmodels`,
  `matplotlib`, `seaborn`, `joblib`, `pyarrow`, `lifelines`.
- **Dashboard**: `streamlit`, `plotly`.
- **Docs**: `mkdocs`, `mkdocs-material`, `pymdown-extensions`, `mkdocstrings[python]`.
- **Dev**: `pytest`, `pytest-cov`, `ruff`.

All pinned in `requirements.txt`, `requirements-dev.txt`, and
`requirements-docs.txt`. The full dependency set also lives in `pyproject.toml`
under `[project]` and `[project.optional-dependencies].dev`.

**Note on real data:** the TCGA-BRCA loader downloads from UCSC Xena on first
use and caches under `data/cache/`, which is gitignored. No authentication is
required. Approximate download size: 40 MB (expression) + 1.7 MB (clinical).

---

## 🤝 Contributing

Contributions are welcome. Please:

1. Fork the repository.
2. Create a feature branch: `git checkout -b feature/my-feature`.
3. Ensure `ruff check src tests` and `pytest -q` pass locally.
4. Commit with a descriptive message.
5. Open a pull request.

The CI workflow will run the full test matrix on your PR. The Docs workflow
will rebuild and redeploy the site when `docs/` or `mkdocs.yml` changes.

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
first — most reviewer questions are answered there. For bugs, please open an
issue on GitHub.

---

## 📄 License

MIT — see [`LICENSE`](LICENSE) for the full text.

---

<p align="center">
  <sub>Built as a demonstration of rigorous, reproducible biomarker-discovery methodology.</sub><br>
  <sub>Synthetic and public TCGA-BRCA data only. Not for clinical use.</sub>
</p>
