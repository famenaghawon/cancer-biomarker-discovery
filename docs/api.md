# API Reference

Public API of the `cancer_biomarkers` package. Everything else is considered
internal and may change between minor versions.

## Configuration

::: cancer_biomarkers.config

### `DataConfig`

Frozen dataclass with all data-source parameters.

| Field | Default | Purpose |
|---|---|---|
| `source` | `"synthetic"` | One of `"synthetic"`, `"tcga_brca"`, `"methylation"` |
| `batch_correct` | `False` | Apply per-batch mean-shift removal before DE |
| `n_samples` | `200` | Number of samples |
| `n_genes` | `500` | Number of genes (or CpG sites) |
| `tumor_ratio` | `0.5` | Fraction of tumour samples |
| `n_true_biomarkers` | `25` | Size of the planted signal panel |
| `effect_size` | `0.5` | Log₂-scale mean shift for true biomarkers |
| `noise_sd` | `1.0` | Per-gene noise multiplier |
| `batch_effect_sd` | `0.25` | Library-size and batch-effect magnitude |
| `random_state` | `42` | Seed for reproducibility |
| `gene_panel` | `()` | Optional custom gene panel (TCGA only) |

### `DEConfig`

| Field | Default | Purpose |
|---|---|---|
| `alpha` | `0.05` | FDR significance threshold |
| `min_abs_log2fc` | `0.6` | Minimum \|log₂FC\| for a significant call |
| `method` | `"fdr_bh"` | Multiple-testing correction method |

### `MLConfig`

| Field | Default | Purpose |
|---|---|---|
| `test_size` | `0.2` | Held-out fraction |
| `n_splits` | `5` | Stratified CV folds |
| `top_k_features` | `30` | Genes selected for ML |
| `random_state` | `42` | Seed for reproducibility |

### `PathConfig`

| Field | Default |
|---|---|
| `results_dir` | `Path("results")` |
| `figures_dir` | `Path("results/figures")` |

`PathConfig.ensure()` creates both directories if missing.

---

## Data sources

### `data_source.synthetic.generate_dataset(cfg) → SyntheticDataset`

Produces a reproducible synthetic expression matrix with known ground-truth
biomarkers, silent tumours, and decoy genes.

### `data_source.tcga.load_tcga_brca(cfg, gene_panel=None) → TCGADataset`

Downloads a curated TCGA-BRCA subset from UCSC Xena, caches it under
`data/cache/`, filters to a gene panel, and returns expression, labels,
clinical metadata, and true biomarkers.

### `multiomics.generate_methylation(cfg) → MethylationDataset`

Generates synthetic methylation β-values with real Illumina 450k probe IDs.

### `multiomics.to_m_values(beta) → DataFrame`

Applies the standard M-value transform:

```
M = log2( (β + ε) / (1 - β + ε) )
```

with ε = 1e-6. Converts β ∈ [0, 1] to an unbounded log-odds scale that the
DE module can process alongside expression data.

---

## Differential expression

### `differential_expression.differential_expression(expression, labels, cfg) → DataFrame`

Per-gene log₂FC, Welch t-test, and BH-FDR. Returns a DataFrame sorted by
q-value with columns:

| Column | Meaning |
|---|---|
| `gene` | Feature name |
| `tumor_mean` | Mean expression in tumour samples |
| `normal_mean` | Mean expression in normal samples |
| `log2fc` | Log₂ fold change |
| `t_stat` | Welch t-statistic |
| `p_value` | Raw two-sided p-value |
| `q_value` | BH-FDR adjusted p-value |
| `significant` | Boolean: `q < alpha` **and** `\|log2fc\| >= min_abs_log2fc` |
| `direction` | `"up"` or `"down"` |

---

## Machine learning

### `ml.train_and_evaluate(expression, labels, cfg, return_test_predictions=False) → dict`

Fast path: one stratified train/test split, one 5-fold CV on the training
features. Returns metrics for both models plus (optionally) test-set
probabilities for ROC plotting.

### `ml.nested_cv_evaluate(expression, labels, cfg) → dict`

Strict path: outer StratifiedKFold for performance, feature selection inside
each fold. Slower, more honest.

### `ml.count_true_in_top_k(top_genes, true_biomarkers) → dict`

Reports how many of the selected top-K genes match the planted biomarker
panel.

---

## Survival

### `survival.cox_ph(expression, time, event, top_genes) → CoxResult`

Fits a Cox proportional hazards model via `lifelines` on standardized
covariates. Automatically drops collinear and zero-variance covariates.

### `survival.kaplan_meier(time, event, stratify_by=None) → list[KMCurve]`

Computes Kaplan–Meier survival curves, optionally stratified by a binary
label.

### `survival.stratify_by_expression(expression, cutoff="median") → Series`

Splits a gene's expression into `"low"` / `"high"` at the given cutoff.

---

## Batch correction

### `batch.correct_batch(expression, batch) → DataFrame`

Per-gene removal of batch mean shifts:

```
corrected[i, g] = expression[i, g] - batch_mean[batch(i), g] + global_mean[g]
```

Location-shift only, not full ComBat. See [Methods](methods.md#batch-effect-correction).

---

## Pipeline

### `pipeline.run(data_cfg=None, de_cfg=None, ml_cfg=None, paths=None) → dict`

End-to-end orchestration. Loads data, optionally corrects batch, runs DE, ML,
survival (if applicable), writes figures, and returns the summary dict that
gets serialised to `results/analysis_summary.json`.
