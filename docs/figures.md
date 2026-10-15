# Figures

This document describes every figure the pipeline emits: what it shows, how it
is built, and how to read it. All figures are written to `results/figures/` by
`src/cancer_biomarkers/plotting.py` during `python scripts/run_pipeline.py`.

All plots use `matplotlib.use("Agg")` so they render headlessly in CI and on
servers without a display.

---

## 1. Volcano plot — `figures/volcano.png`

### What it shows

Every gene plotted as a point, with:

- **x-axis** — log₂ fold change (tumour vs. normal)
- **y-axis** — −log₁₀(q-value), i.e. statistical significance (higher = more significant)

### How it is built

```python
df["-log10q"] = -np.log10(df["q_value"].clip(lower=1e-300))
sig = df["significant"]
```

Three scatter groups:

| Group | Colour | Condition |
|---|---|---|
| Non-significant | Grey | `not significant` |
| Upregulated | Crimson | `significant and log2fc > 0` |
| Downregulated | Steel blue | `significant and log2fc < 0` |

Threshold lines:

- Horizontal dashed line at `y = -log10(alpha)` (default `alpha = 0.05` → y ≈ 1.3)
- Vertical dashed lines at `x = ±min_abs_log2fc` (default `0.6`)

The title embeds both thresholds explicitly:

```
Volcano plot — tumor vs normal (|log2FC| >= 0.6, q < 0.05)
```

### How to read it

- Points **inside** the central cross (|x| < `min_abs_log2fc` **or** y < `-log10(alpha)`)
  are non-significant under the DE rule. They should be grey.
- Points in the **upper-left** quadrant are significantly downregulated.
- Points in the **upper-right** quadrant are significantly upregulated.
- The horizontal spread of a point's y-value shows how much statistical evidence
  supports the effect; a small but consistent log₂FC can be significant.

### How to sanity-check the image

1. Count the coloured (red + blue) points. This count should equal
   `n_significant_biomarkers` in `analysis_summary.json`.
2. The leftmost coloured point's x-value should be ≈ −max(|log₂FC|) among
   significant genes; the rightmost should be ≈ +max(|log₂FC|).
3. The horizontal dashed line should sit at y ≈ 1.3 for α = 0.05.
4. Verify the true biomarker panel (real gene symbols) appears predominantly
   in the coloured region.

### Known edge cases

- **Zero significant genes** — plot is entirely grey. Not an error; it means
  the effect size is too small for the chosen thresholds. Reduce
  `min_abs_log2fc` or raise `DataConfig.effect_size`.
- **Heavy tail** — if one or two genes have very high |log₂FC|, the x-range
  stretches and the central cluster collapses. This is expected; the y-axis
  is still informative.

---

## 2. Heatmap — `figures/heatmap.png`

### What it shows

Top-25 biomarkers (by |Welch t-statistic| on the training fold) as rows, samples
as columns, with expression values z-scored per gene.

### How it is built

```python
sub = expression[genes].copy()
if log_transform:
    sub = np.log2(sub.clip(lower=1e-9))
mat = sub.T
mat = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, 1), axis=0)
order = labels.sort_values(kind="stable").index
mat = mat.loc[:, order]
```

Two transformations are applied:

1. **log₂-transform** — compresses the dynamic range so a few high-expression
   genes do not dominate the colour scale. `clip(lower=1e-9)` guards against
   log of zero.
2. **z-score per gene** — subtract mean, divide by std, so each row is on the
   same scale. Without this, a gene with expression ~8000 would swamp a gene
   with expression ~80, and both patterns would be invisible.

Samples are ordered by class (`labels.sort_values`) so tumour and normal samples
cluster in two blocks.

### Colour scheme

- **Diverging colourmap (`vlag`)** with `center = 0.0`
- **Red** — above-average expression for that gene
- **Blue** — below-average expression for that gene
- **White** — at the gene's mean

The colorbar is labelled "z-score (per gene)" so the reader knows the units are
not raw expression.

### How to read it

- A **horizontal stripe** (consistent colour across columns) means the gene is
  uniformly high or low across samples — no biological contrast.
- A **checkerboard split** — red on the tumour half, blue on the normal half,
  or vice versa — means the gene separates the two classes. These are your
  candidates.
- A **noisy row** with no pattern means the gene ranks in the top-25 by
  t-statistic but the effect is driven by outliers; worth inspecting the
  individual values in `de_results.csv`.

### How to sanity-check the image

1. The leftmost block of columns (normals, or tumours depending on label
   ordering) should show roughly opposite colour patterns to the rightmost
   block for the true biomarkers.
2. The row labels should be recognisable real gene symbols (BRCA1, TP53, EGFR,
   ...) plus possibly some `Gene_XXXX` fillers, since the top-25 by t-statistic
   may include decoys.
3. Number of rows = 25 (or fewer if fewer genes were selected; `top_k_features`
   default is 30, but the figure takes `[:25]`).

### Known edge cases

- **Uniformly coloured rows** — the gene was selected but the effect is subtle.
  This is honest; not every top-25 gene will show a clean split.
- **Diagonal gradient** — indicates a batch or library-size effect that the
  pipeline does not correct. Our generator injects both, so a mild gradient is
  expected. It does not invalidate the class separation.

---

## 3. ROC curve — `figures/roc_curve.png`

### What it shows

Receiver-operating-characteristic curve for logistic regression on the
held-out test set.

- **x-axis** — false positive rate (1 − specificity)
- **y-axis** — true positive rate (sensitivity)
- **Diagonal dashed line** — random classifier (AUC = 0.5)
- **Solid curve** — the fitted classifier
- **Legend label** — `logistic_regression (AUC = <value>)`

### How it is built

```python
X_tr, X_te, y_tr, y_te = train_test_split(..., stratify=labels, random_state=seed)
genes = _select_top_k(X_tr, y_tr, top_k_features, seed)
model = LogisticRegression(C=0.1, max_iter=2000).fit(X_tr[genes], y_tr)
y_proba = model.predict_proba(X_te[genes])[:, 1]
fpr, tpr, _ = roc_curve(y_te, y_proba)
auc = roc_auc_score(y_te, y_proba)
```

The same split and same feature selection as the reported ML metrics, so the
AUC in the legend matches `ml.models.logistic_regression.roc_auc` in
`analysis_summary.json` exactly.

### How to read it

- AUC near 1.0 → perfect separation (suspicious for a biomarker task)
- AUC 0.80–0.95 → strong signal, realistic
- AUC 0.60–0.80 → weak signal, still useful
- AUC near 0.5 → classifier is at chance

The shape matters too:

- A curve that rises steeply at the left (low FPR) and flattens near the top
  means the classifier confidently identifies tumours without many false
  positives — good.
- A curve that follows the diagonal for the first half then jumps means the
  classifier is unsure at low thresholds — the AUC will still look okay but
  the operating point choice is critical.

### How to sanity-check the image

1. The AUC printed in the legend must match the JSON value for
   `logistic_regression.roc_auc` (rounding aside).
2. The curve must be monotonically increasing.
3. The endpoints must be (0, 0) and (1, 1) — any ROC implementation guarantees
   this.

### Known edge cases

- **AUC = 1.0** — the classifier separates perfectly. On our synthetic data
  with `effect_size = 0.5` this should not happen; if it does, lower
  `DataConfig.effect_size` to 0.4 and re-run.
- **AUC < 0.6** — the effect is too subtle; raise `effect_size` or reduce the
  silent-tumour fraction in `data.py`.

---

## 4. Reproducing the figures

Figures are not committed to the repository (see `.gitignore`). To regenerate
them:

```bash
python scripts/run_pipeline.py
```

Then inspect `results/figures/`. Each figure is 150 DPI PNG, roughly 40–80 KB.

To view a specific figure from the command line on Windows:

```powershell
Start-Process results\figures\volcano.png
Start-Process results\figures\heatmap.png
Start-Process results\figures\roc_curve.png
```

To regenerate with different thresholds without editing code:

```python
from cancer_biomarkers.config import DataConfig, DEConfig, MLConfig
from cancer_biomarkers.pipeline import run

run(
    data_cfg=DataConfig(effect_size=0.7),   # stronger signal
    de_cfg=DEConfig(min_abs_log2fc=0.8),    # stricter DE
    ml_cfg=MLConfig(top_k_features=50),     # more features
)
```

The plots will reflect the new thresholds automatically — no figure code changes
required.

---

## 5. What is *not* plotted

The following are deliberately absent, to avoid implying analytical steps the
pipeline does not perform:

- **Survival curves (Kaplan–Meier).** The clinical file contains a `stage`
  column, but no survival times. Plotting KM curves without time-to-event data
  would be misleading.
- **PCA / UMAP / t-SNE.** Not needed; the pipeline's message is DE + ML, not
  dimensionality reduction. If you want a UMAP of the expression matrix as a
  sanity check, that is a valid addition — but it is not part of the pipeline's
  scientific claim.
- **Confusion matrix heatmap.** The confusion matrix is stored in
  `analysis_summary.json` under `ml.models.<name>.confusion_matrix`; it is not
  plotted because the same information is captured by the accuracy, precision,
  and recall metrics in the summary.
- **Feature importance plot.** Random forest importances are not currently
  exported. Adding one is straightforward and would complement the heatmap.

If you want any of these, each is a ~15-line addition to `plotting.py` and a
one-line call in `pipeline.py`. Prioritise the survival curves only after
adding real time-to-event data to the generator.
