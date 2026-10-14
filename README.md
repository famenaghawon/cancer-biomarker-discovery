# 🧬 Cancer Biomarker Discovery Pipeline

> Educational, end-to-end pipeline for cancer biomarker discovery on **synthetic TCGA-like** data.
> Combines rigorous differential expression (Welch t-test + BH-FDR) with a leakage-safe
> machine-learning workflow.

## ⚠️ Disclaimer
This is a **demo / research tool**. All data is synthetic and results are **not for clinical use**.

## ✨ What it does

1. Generates a reproducible synthetic expression matrix with known ground-truth biomarkers.
2. Runs differential expression: log₂ fold change, Welch t-test, Benjamini–Hochberg FDR.
3. Trains classification models with **feature selection on training data only** (no leakage).
4. Reports cross-validated ROC-AUC in addition to a held-out test set.
5. Emits volcano plot and heatmap, plus a JSON summary and CSV biomarker table.

## 🚀 Quick Start


git clone https://github.com/famenaghawon/cancer-biomarker-discovery.git
cd cancer-biomarker-discovery

python -m venv .venv
# Windows PowerShell
.venv\Scripts\Activate.ps1
# macOS/Linux
source .venv/bin/activate

pip install -r requirements.txt

python scripts/run_pipeline.py

Results land in results/:


results/
├── analysis_summary.json
├── clinical_data.csv
├── de_results.csv
├── significant_biomarkers.csv
├── expression_full.parquet
├── ml_results.pkl
└── figures/
    ├── volcano.png
    └── heatmap.png

🧪 Tests

pip install -r requirements-dev.txt
pytest


📁 Project Structure

cancer_analysis/
├── src/cancer_biomarkers/
│   ├── __init__.py
│   ├── config.py
│   ├── data.py
│   ├── differential_expression.py
│   ├── logging_utils.py
│   ├── ml.py
│   ├── pipeline.py
│   └── plotting.py
├── scripts/run_pipeline.py
├── tests/
├── results/
├── requirements.txt
├── requirements-dev.txt
├── pyproject.toml
├── README.md
├── PROJECT_SUMMARY.md
└── LICENSE


🔬 Methods

Differential expression

· Per-gene log₂FC = log₂((tumor_mean + ε) / (normal_mean + ε))
· Welch's t-test (unequal variances)
· Benjamini–Hochberg FDR
· Significant: q < 0.05 and |log₂FC| ≥ 1

Machine learning

· Models: logistic regression, random forest
· Feature selection is done on the training fold only (avoids leakage)
· 5-fold stratified cross-validation on the training set
· Held-out test set for final metrics

Synthetic data

· Baseline expression ~ N(8, 1.5)
· Per-gene Gaussian noise + per-sample batch effect
· Tumor effect on a small set of "true" biomarkers (log₂-scale shift)
· 30 real gene symbols (BRCA1, TP53, EGFR, …) reused as ground-truth markers

📊 Reported Metrics

Metrics are written to results/analysis_summary.json and reflect the actual run.
They will vary slightly with DataConfig/MLConfig. The pipeline deliberately uses a
small effect size so models do not trivially hit 100%.

📝 Citation


@software{cancer_biomarker_2025,
  title  = {Cancer Biomarker Discovery Pipeline},
  author = {Famenaghawon, Freedom},
  year   = {2025},
  url    = {https://github.com/famenaghawon/cancer-biomarker-discovery}
}


📧 Contact

f.e.amenaghawonfreedom@gmail.com

📄 License

MIT — see LICENSE.
