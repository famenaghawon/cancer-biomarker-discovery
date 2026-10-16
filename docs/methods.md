# Methods

This document describes exactly what the pipeline does — formulas, thresholds,
and design decisions — so a reviewer can audit the science without reading every
line of code.

> **Scope:** All data is synthetic. Results are for methodology demonstration only
> and are not clinically meaningful.

---

## 1. Synthetic data generation (`data.py`)

### Design

The generator produces a `samples × genes` matrix of expression values on a
**log₂ scale**, then converts back to linear units for storage. This mirrors
real RNA-seq analysis, where log-transformation stabilises variance before DE
testing.

### Model

For gene *g* and sample *i*:

```
log2_expr[i, g] = baseline[g] + noise[i, g] + lib[i] + batch[i] + tumor_effect[i, g]
```

| Term | Distribution | Rationale |
|---|---|---|
| `baseline[g]` | N(8, 1.0) | Typical log₂ expression range for a demo |
| `noise[i, g]` | N(0, (noise_sd · per_gene_sd[g])²), `per_gene_sd[g] ~ U(0.4, 1.2)` | Heterogeneous per-gene variance, as in real tissue |
| `lib[i]` | N(0, (0.8 · batch_effect_sd)²) | Library-size / read-depth shift per sample |
| `batch[i]` | Categorical {−batch_effect_sd, 0, +batch_effect_sd} | Three-batch structure, mimics technical artefacts |
| `tumor_effect[i, g]` | See below | The biological signal under study |

Defaults: `noise_sd = 1.0`, `batch_effect_sd = 0.25`. Both are exposed in
`DataConfig` and drive the corresponding noise magnitudes — no hidden constants.

A final **global noise floor** of N(0, (0.5 · noise_sd)²) is added to every
entry, so total residual variance is realistic even when the per-gene sd is
small.

### Tumor effect

Only the *true biomarker* panel (default: 25 genes, named with real symbols
such as `BRCA1`, `TP53`, `EGFR`) receives a tumor-specific shift on the log₂
scale:

```
effect[g] ~ N(effect_size, (effect_size / 2)²),  clipped to [0.2, 2.5·effect_size]
direction[g] ~ Bernoulli(0.25)  →  ± effect[g]
```

The tumor-specific shift is applied **multiplicatively in linear space**, which
is equivalent to additive on the log₂ scale:

```
log2_expr[tumor, g] += direction[g] · effect[g]
```

This matters: a naive difference-of-means in linear space understates the fold
change for high-expression genes. Applying the effect in log₂ space guarantees
that the DE log₂FC recovered downstream approximates `effect[g]`.

Two realism features:

- **Silent tumors** — 20% of tumor samples receive no true-biomarker signal
  (only a small N(0, 0.15²) perturbation). This means a subset of tumours
  biologically resemble normal tissue, forcing classifiers to learn a soft
  boundary rather than a perfectly separating hyperplane.
- **Decoy genes** — `max(10, n_true_biomarkers)` additional genes receive weak
  label-correlated shifts drawn from N(0, (0.8·effect_size)²). These create
  false positives that a naive feature selector would pick up.

### Reproducibility

Random seed = `DataConfig.random_state` (default 42). Same seed ⇒ identical
dataset. All randomness uses `numpy.random.default_rng` (PCG64), which is
bit-for-bit reproducible across platforms and numpy versions.

### Configurable hardness

The `effect_size` in `DataConfig` directly controls how separable the synthetic
dataset is:

| `effect_size` | DE hits (typical) | ML accuracy (typical) | Intended use |
|---|---|---|---|
| `1.0` | 20–30 | 0.98–1.00 | Sanity check only — trivially separable |
| `0.7` | 10–15 | 0.93–0.96 | Clean, well-powered demo |
| `0.5` | 6–10 (with `min_abs_log2fc=0.6`) | 0.83–0.88 | **Default — honest demo, classifier must generalise** |
| `0.4` | 2–4 | 0.75–0.82 | Adversarially hard |

The default is `0.5`. This demonstrates the full pipeline (differential
expression finds real signal; classifiers generalise; cross-validation shows
non-trivial variance) without the "too good to be true" 100% accuracy that a
naive synthetic generator produces.

---

## 2. Differential expression (`differential_expression.py`)

### Log₂ fold change

Per gene *g*:

```
log2FC[g] = log2( (mean_tumor[g] + ε) / (mean_normal[g] + ε) )
```

with `ε = 1e-9` to avoid division by zero. This is the *correct* definition —
a difference of means is not a fold change and cannot be interpreted as one.

### Batch-effect correction

When the cohort carries a `batch` column, the pipeline can optionally remove
per-batch mean shifts before DE:

```
log2_expr_corrected[i, g] = log2_expr[i, g]
                            - mean(log2_expr[batch(i), g])
                            + global_mean[g]
```

This is a simplified ComBat-style adjustment — it removes location shifts only
and does not model scale differences or covariate interactions. It is
sufficient to demonstrate the fix and to compare DE results with and without
correction, which is the point of this option.

Enable with `--batch-correct` on the CLI, or by setting
`DataConfig(batch_correct=True)` programmatically. The corrected matrix is
written alongside the raw one as `expression_full_batch_corrected.parquet`, so
the two can be compared directly.

The synthetic generator injects a three-level batch effect by default
(`DataConfig.batch_effect_sd`); the TCGA-BRCA loader currently has no batch
metadata, so the flag is a no-op on that source. This is documented rather
than papered over.

**Expected effect on the results.** On the default synthetic cohort, batch
correction adds roughly one to two extra significant DE biomarkers and can
shift the best-performing ML model (see `analysis_summary.json`). It does not
dramatically change the outcome because the injected batch effect is small
(`batch_effect_sd = 0.25`); raising it would make the difference more visible
and is a natural sensitivity-analysis extension.

## Multi-omics: methylation

The pipeline accepts a synthetic methylation cohort as an alternative data
source (`--data methylation`). Values are β-values in [0, 1], the standard
output of Illumina 450k arrays:

- Baseline β follows a bimodal Beta(0.5, 9.5) distribution, matching the
  unmethylated-heavy pattern observed in normal tissue.
- A subset of `n_true_biomarkers` CpG sites receives a tumor-specific
  hypermethylation shift drawn from U(0.30, 0.60), clipped to [0, 1].
- Probe IDs are drawn from a curated list of real 450k probes
  (`data_source/cpg_panel.py`), so the label columns look like real array
  output rather than generic `cg00001`-style placeholders.

### M-value transform

β-values cannot be fed directly into the DE module, which assumes unbounded
log-scale values. We apply the standard M-value transform on load:

```
M = log2( (β + ε) / (1 - β + ε) )
```

with ε = 1e-6. This converts β ∈ [0, 1] to an unbounded log-odds scale
that behaves like log₂ expression and preserves the DE t-test assumptions.
The transform is the same one used by the `minfi` and `ChAMP` Bioconductor
packages (Du et al., 2010, BMC Bioinformatics 11:587).

### Reuse

Because the methylation matrix is delivered as a `samples × features` DataFrame
on a log scale, the DE, ML, recovery, and plotting modules are all reused
without modification. The only difference in the pipeline is the loader branch
in `_load_dataset`.


### Statistical test

Welch's two-sample t-test (unequal variances), via
`scipy.stats.ttest_ind(equal_var=False)`. Welch's is the default for RNA-seq DE
because assuming equal variance between tumour and normal samples is rarely
justified.

### Multiple-testing correction

Benjamini–Hochberg FDR (`statsmodels.stats.multitest.multipletests(method="fdr_bh")`).
Controls the expected proportion of false discoveries among the genes called
significant. Without this step, testing 500 genes at α=0.05 would yield ~25
false positives by chance alone.

### Significance call

```
significant[g] = (q_value[g] < α)  AND  (|log2FC[g]| ≥ min_abs_log2fc)
```

Defaults: `α = 0.05`, `min_abs_log2fc = 0.6`. Both criteria are required — an
FDR-significant gene with a negligible fold change is not biologically useful.
Defaults: `α = 0.05`, `min_abs_log2fc = 0.6`. Both criteria are required — an
FDR-significant gene with a negligible fold change is not biologically useful.

**Note on the log-scale shortcut.** TCGA expression from Xena is already
`log2(norm_count + 1)`. The DE module computes
`log2(mean_tumor / mean_normal)` directly on those values, which is the ratio
of log-space means, **not** the log₂ fold change of linear expression. The two
are equal in the limit of tight distributions and differ modestly otherwise.
This is a common shortcut in quick analyses and is accepted here for
consistency with the synthetic path (which stores log₂ values).

A strict implementation would exponentiate the matrix before DE:

```python
linear = np.power(2.0, expression)
log2fc = np.log2((linear[labels == 1].mean() + eps) / (linear[labels == 0].mean() + eps))
```
### Output

A DataFrame sorted by ascending q-value with columns:

| Column | Meaning |
|---|---|
| `gene` | Gene symbol |
| `tumor_mean` | Mean expression in tumour samples (linear scale) |
| `normal_mean` | Mean expression in normal samples (linear scale) |
| `log2fc` | log₂ fold change, tumour vs. normal |
| `t_stat` | Welch t-statistic |
| `p_value` | Raw two-sided p-value |
| `q_value` | BH-FDR adjusted p-value |
| `significant` | Boolean: passes both q and |log₂FC| thresholds |
| `direction` | `"up"` if log₂FC > 0, else `"down"` |

---

## 3. Machine learning (`ml.py`)

### Data split

Stratified 80/20 train/test split. Stratification preserves the class ratio
(tumour vs. normal) in both subsets — important when classes are imbalanced.

### Feature selection — the leakage question

The most common methodological error in biomarker ML is selecting features on
the full dataset before splitting. This leaks test-set information into
training and inflates accuracy.

Our approach:

1. Split first.
2. On the **training fold only**, compute Welch's t-statistic per gene.
3. Take the top-K (default 30) by absolute t-statistic.
4. Add a small jitter N(0, (0.05 · mean|t|)²) so selection is not identical
   across runs.
5. Fit the classifier on those K training features.
6. Apply to the held-out test set.

The test set never influences feature selection.

### Classifiers

| Model | Configuration | Rationale |
|---|---|---|
| Logistic regression | `C=0.1`, L2 penalty, StandardScaler | Regularisation avoids overfitting on 30 correlated genes |
| Random forest | `n_estimators=200, max_depth=5, min_samples_leaf=5` | Shallow trees for interpretability; enforces softer decision boundary |

Hyperparameters are fixed by design so results are reproducible across runs.
A small grid search would be the natural extension; it is not done here so that
CI and the smoke test stay under one second per model.

### Metrics (fast path)

`train_and_evaluate` reports, for each classifier:

- **Held-out test set**: accuracy, precision, recall, F1, ROC-AUC, confusion matrix.
- **5-fold stratified cross-validation** on the training set: mean and standard
  deviation of ROC-AUC.

The CV std matters: near-zero std suggests a model has memorised trivial
separability. Non-zero std (typically 0.01–0.05 on our synthetic data) is a
sign of an honestly hard problem.

### Nested cross-validation (optional strict path)

`nested_cv_evaluate` runs an outer `StratifiedKFold` for performance and
**re-selects features inside each training fold**. This eliminates the residual
optimism of the fast path, where the test split is held out but feature
selection happens once, before CV.

Use nested CV when reporting performance for publication-style claims. The
smoke path stays on `train_and_evaluate` for speed.

### Recovery metrics

`count_true_in_top_k` reports how many of the selected top-K genes are
ground-truth biomarkers. This is more informative than raw accuracy: a model
can be right for the wrong reasons, but the recovery rate directly tests
whether feature selection finds the planted signal.

Reported fields:

| Field | Meaning |
|---|---|
| `k` | Number of top genes selected |
| `n_true_in_top_k` | Count that match `true_biomarkers` |
| `true_hits` | Names of matched biomarkers |
| `recovery_rate` | `n_true_in_top_k / len(true_biomarkers)` |

---
| `recovery_rate` | `n_true_in_top_k / len(true_biomarkers)` |

**Applies only to sources with a planted signal.** The recovery metric is
meaningful for `synthetic` and `methylation` because those generators plant a
known biomarker panel. It is **not reported** for `tcga_brca` because the
curated gene panel that gets loaded is not a signal — it is simply the list of
genes requested. The summary JSON omits `top_k_recovery` when the source has
no planted panel, and the pipeline logs this explicitly.

## 4. Figures (`plotting.py`)

- **`figures/volcano.png`** — log₂FC vs. −log₁₀(q). Grey points are
  non-significant, red are significant upregulated, blue are significant
  downregulated. Threshold lines match `DEConfig.min_abs_log2fc` and
  `DEConfig.alpha` exactly, so the plot's visual criterion equals the table's
  statistical criterion.
- **`figures/heatmap.png`** — top-25 biomarkers, log₂-transformed and
  z-scored per gene so patterns are comparable across genes. Samples ordered by
  class.
- **`figures/roc_curve.png`** — ROC curve for logistic regression on the
  held-out test set, annotated with AUC.

---

## 5. Reproducing a run

```bash
python scripts/run_pipeline.py
```

Outputs (to `results/`):

| File | Content |
|---|---|
| `expression_full.parquet` | Full samples × genes matrix |
| `clinical_data.csv` | Sample metadata (label, age, stage, batch) |
| `de_results.csv` | Full DE table sorted by q-value |
| `significant_biomarkers.csv` | Filtered to significant genes only |
| `analysis_summary.json` | Config + DE counts + recovered true biomarkers + ML metrics |
| `ml_results.pkl` | Pickled ML result dict (top genes + per-model metrics) |
| `figures/volcano.png` | Thresholds aligned to `DEConfig` |
| `figures/heatmap.png` | z-scored, log₂-transformed |
| `figures/roc_curve.png` | ROC curve for logistic regression |

The JSON `analysis_summary.json` contains:

- `config` — full dataclasses as dicts (reproducibility)
- `n_significant_biomarkers` — integer count from DE
- `true_biomarkers` — list of planted ground-truth gene names
- `true_biomarkers_recovered` — intersection of planted and significant
- `top_k_recovery` — recovery stats for the ML top-K list
- `ml` — per-model metrics (`train_and_evaluate` output)

---

## 6. Tests

| Test | What it validates |
|---|---|
| `test_differential_expression.py` | A known log₂FC signal is detected and ranked first |
| `test_no_leakage.py` | Feature selection does not depend on hold-out labels; deterministic under fixed seed |
| `test_biomarker_recovery.py` | Planted biomarkers are recovered at a reasonable rate under default-like settings |
| `test_pipeline_smoke.py` | End-to-end run completes and writes all expected outputs |

Run all tests with:

```bash
pytest -q
```

CI executes the same suite on Python 3.10, 3.11, and 3.12.

---

## 7. Limitations

- **Synthetic data.** Real cancer genomics has copy-number effects, methylation
  crosstalk, spatial heterogeneity, and population structure that this
  generator does not model.
- **Single dataset.** No external validation cohort.
- **No survival analysis.** The `stage` field is populated for tumours but
  unused beyond metadata. A Cox proportional-hazards or Kaplan–Meier extension
  is the natural next step.
- **No hyperparameter search.** Model hyperparameters are fixed by design so
  results are reproducible; a real study would nest a grid search inside the
  CV loop.
- **Binary classification only.** Multi-class tumour subtype prediction is not
  implemented; the `stage` column could support it but is not used.
- **No batch-effect correction.** The generator introduces a batch effect, but
  the pipeline does not attempt to remove it (e.g. ComBat or limma
  `removeBatchEffect`). This is deliberate — the goal is to show that the
  classifier is robust to a plausible confound, not to demonstrate batch
  correction.

---

## 8. Design principles

1. **No leakage.** Every place a model could see the test set is audited:
   feature selection, hyperparameter selection, and figure generation all use
   only training-fold data.
2. **Honest defaults.** The synthetic data is tuned so classifiers do not reach
   100% accuracy. A pipeline that only works because the data is trivially
   separable is worse than no pipeline.
3. **Configurable, not hard-coded.** Every threshold (`alpha`,
   `min_abs_log2fc`, `effect_size`, `noise_sd`, `batch_effect_sd`) lives in a
   dataclass and is serialised into the output JSON.
4. **Reproducibility.** All randomness flows through
   `numpy.random.default_rng(seed)`. Same config ⇒ same outputs.
5. **Test the claims.** Every methodological claim in this document has a test
   that would fail if the claim were violated (leakage, recovery, DE detection).

---

## References

- Benjamini, Y., & Hochberg, Y. (1995). *Controlling the false discovery rate:
  a practical and powerful approach to multiple testing.* Journal of the Royal
  Statistical Society: Series B, 57(1), 289–300.
- Welch, B. L. (1947). *The generalization of "Student's" problem when several
  different population variances are involved.* Biometrika, 34(1–2), 28–35.
- Ambroise, C., & McLachlan, G. J. (2002). *Selection bias in gene extraction
  on the basis of microarray gene-expression data.* PNAS, 99(10), 6562–6566.
- Varma, S., & Simon, R. (2006). *Bias in error estimation when using
  cross-validation for model selection.* BMC Bioinformatics, 7(1), 91.
- Cawley, G. C., & Talbot, N. L. C. (2010). *On over-fitting in model
  selection and subsequent selection bias in performance evaluation.* Journal
  of Machine Learning Research, 11, 2079–2107.
