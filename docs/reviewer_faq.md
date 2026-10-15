# Reviewer FAQ

Answers to the questions a skeptical reviewer, hiring manager, or peer
researcher is likely to ask when reading this repository. Each answer points at
the concrete file or line of code that supports the claim, so nothing has to be
taken on trust.

---

## Q1. Why synthetic data? Why not real TCGA data?

**Short answer:** network access, licensing, and reproducibility.

**Long answer:** Real TCGA data is publicly available but:

- Requires an authenticated download (GDC API keys, or a mirror with its own
  terms). CI cannot fetch it reliably, and the repo would be unusable offline.
- Has licence terms that vary by cohort and downstream use; embedding a real
  subset in a demo repo risks violating redistribution terms.
- Is large enough that a fresh clone plus `pytest` would take minutes instead
  of seconds.

The synthetic generator models the features that matter for a DE + ML pipeline:

- log₂-scale expression with per-gene variance (`data.py`)
- batch and library-size effects
- heterogeneous effect sizes and direction on a known biomarker panel
- silent tumours (a fraction of tumours that look normal on the panel)
- decoy genes that create plausible false positives

This gives the pipeline a *known ground truth*, which real data cannot. The
recovery rate (`top_k_recovery` in `analysis_summary.json`) directly measures
whether the pipeline found the planted signal, not just whether the classifier
happened to fit.

Swapping in real data is the subject of a documented future extension in
`docs/methods.md` §7.

---

## Q2. Isn't the classifier still suspiciously accurate at ~85%?

**Short answer:** Yes, and that is the point. It is not 100%.

**Long answer:** With `effect_size = 0.5` on the log₂ scale and 20% silent
tumours, an honest classifier should sit in the 0.80–0.90 accuracy range. If
you see 0.95 or above, the synthetic data has drifted toward trivial
separability — three knobs restore realism:

1. `DataConfig.effect_size` — lower to `0.4` (`config.py`)
2. `data.py` — raise `per_gene_sd` from `U(0.4, 1.2)` to `U(0.8, 2.0)`
3. `data.py` — raise the silent-tumour fraction from `0.20` to `0.30`

The hardness table in `docs/methods.md` §1 documents the expected range for
each `effect_size` so you can reproduce it in seconds.

Near-perfect accuracy is a *warning sign* in biomarker ML, not a goal. A
classifier that separates perfectly on synthetic data has learned something
about the generator, not about biology.

---

## Q3. How do I know there is no data leakage?

**Short answer:** Two code paths are audited and one test enforces it.

**Long answer:**

- **Feature selection** is computed on `X_train` only. See
  `src/cancer_biomarkers/ml.py` — `_select_top_k(X_train, y_train, ...)`.
  The test fold is never passed to `_select_top_k`.
- **Cross-validation** is run on the *already selected* training features
  (fast path) or re-selects inside each fold (nested path,
  `nested_cv_evaluate`). Neither path touches the test set during selection.
- **Test-set metrics** are computed once at the end on `X_test_k`, which
  was never fed to any fitting or selection step.
- **Test enforcement** — `tests/test_no_leakage.py::test_select_top_k_ignores_holdout_labels`
  constructs a dataset, mutates the hold-out labels, and asserts that the
  selected top-K does not change. If selection were contaminated by hold-out
  information, this test would fail.

The nested CV path exists because the fast path still has a small optimistic
bias from selecting once before the outer CV loop. The methods doc
(`docs/methods.md` §3) explains when to use which.

---

## Q4. Why Welch's t-test and not a negative-binomial GLM (as DESeq2 or edgeR do)?

**Short answer:** The generator works on log₂-scale data, not raw counts.

**Long answer:** DESeq2 and edgeR model **raw count data** with a
negative-binomial distribution to account for mean-variance dependence. Our
generator emits continuous log₂-scale values, converted to linear for storage.
For continuous data, Welch's t-test is the standard test — it does not assume
equal variances between groups, which is usually a bad assumption in tumour vs.
normal comparisons.

If you swap in real count data, the correct move is to fit DESeq2 or edgeR
directly. That would require a different data load path (counts, not
continuous values) and a different DE module. It is a genuine extension, not a
drop-in replacement.

---

## Q5. Why Benjamini–Hochberg, and not Bonferroni?

**Short answer:** BH controls the false discovery *rate*; Bonferroni controls
the family-wise *error rate*. For genomics, BH is the field standard.

**Long answer:** Bonferroni divides α by the number of tests, which is
extremely conservative when you have hundreds or thousands of genes. You would
miss real signal. BH controls the expected *proportion* of false discoveries
among rejected hypotheses, which is what you actually care about in a
discovery setting: "of the genes I call significant, at most X% are expected to
be false positives."

BH is implemented via `statsmodels.stats.multitest.multipletests(method="fdr_bh")`
in `differential_expression.py`. It is the same correction used by DESeq2,
edgeR, and limma for standard RNA-seq DE.

If you want a more conservative test, `multipletests` also supports
`method="bonferroni"` — one-line change, but not the default for good reason.

---

## Q6. Why log₂ fold change and not just "tumour mean / normal mean"?

**Short answer:** log₂FC is symmetric, interpretable as "doubling / halving,"
and is the field convention.

**Long answer:**

- A ratio of 2 in one direction and 0.5 in the other are symmetric in log₂
  (+1 and −1), but asymmetric in linear scale (2 and 0.5). Symmetry matters
  when you average across genes or rank by effect size.
- "Doubling of expression" maps to a clean log₂FC = 1. That is why the
  default `min_abs_log2fc = 0.6` corresponds to roughly a 1.5× change.
- A **difference of means** (which the original version of this repo used) is
  not a fold change at all. Two genes with the same difference-of-means can
  have wildly different fold changes depending on baseline expression. On a
  log₂ scale, differences correspond to fold changes; on a linear scale, they
  do not.

See `differential_expression.py`:
`log2fc = np.log2((tumor_mean + _EPS) / (normal_mean + _EPS))`.

---

## Q7. How is the "top-K recovery" metric computed, and why does it matter?

**Short answer:** It counts how many of the top-K selected genes are in the
planted true-biomarker panel.

**Long answer:** Accuracy tells you whether the classifier *works*, but not
*whether it found the right thing*. A model can be 95% accurate while looking
at the wrong genes. Recovery directly answers: "of the 30 genes the pipeline
selected, how many are the ones we planted?"

Implementation (`ml.py::count_true_in_top_k`):

```python
true_set = set(true_biomarkers)
hits = [g for g in top_genes if g in true_set]
return {
    "k": len(top_genes),
    "n_true_in_top_k": len(hits),
    "true_hits": hits,
    "recovery_rate": len(hits) / max(1, len(true_set)),
}
```

The result is written to `analysis_summary.json` under `top_k_recovery` and
logged at the end of the pipeline. A companion test
(`tests/test_biomarker_recovery.py`) asserts that at least half of the true
biomarkers appear among the significant DE genes under default-like settings —
if feature selection were broken, this test would fail.

---

## Q8. Why fixed hyperparameters instead of a grid search?

**Short answer:** Reproducibility and CI speed.

**Long answer:** A grid search adds a nested loop that changes the selected
hyperparameters based on the fold split. Two runs of the same script could
produce different models, which makes the pipeline's outputs
non-deterministic from a reviewer's perspective. It also multiplies the run
time by the number of grid points (often 5–20×), and the CI runs three Python
versions.

The current design uses:

- Logistic regression `C = 0.1` (mild L2 regularisation; enough for 30
  features without overfitting)
- Random forest `max_depth = 5, min_samples_leaf = 5, n_estimators = 200`
  (shallow trees, no overfitting)

`nested_cv_evaluate` in `ml.py` is the natural home for a hyperparameter
search — the inner loop of nested CV is exactly where the search should live.
Adding `GridSearchCV` there is a ~10-line change and would not affect the
smoke-path CI.

The trade-off is documented in `docs/methods.md` §7 (Limitations).

---

## Q9. Why two classifiers, and not a deep learning model?

**Short answer:** The right model for the data, and the right model for the
reader.

**Long answer:**

- **Data size.** 200 samples × 500 genes is well within the regime where
  regularised linear models (logistic regression) and shallow ensembles
  (random forest) are competitive with deep models. A neural network on 200
  samples would overfit without heavy regularisation and extensive tuning.
- **Interpretability.** Logistic regression coefficients directly indicate
  which genes push toward tumour vs. normal. Random forest importances give a
  complementary, non-linear view. Both are straightforward to explain in a
  methods section. A deep model's "importance" is a post-hoc approximation.
- **Domain convention.** The biomarker-discovery literature overwhelmingly
  uses linear / tree-based models for expression-based classification. A demo
  that mirrors that convention is more useful than one that chases a trend.

Adding a third model (e.g. SVM, gradient boosting) is a two-line change in
`ml._make_models`. It is not added by default because two classifiers already
demonstrate the pipeline's modelling surface and keep CI fast.

---

## Q10. What would you do differently if this were a real study?

**Short answer:** Six concrete changes, in priority order.

**Long answer:**

1. **Real data with external validation.** Download a well-characterised cohort
   (e.g. TCGA-BRCA), split by *institution* or *batch*, not randomly. Random
   splits leak technical information across folds.
2. **Batch-effect correction.** Run ComBat or limma's `removeBatchEffect`
   before DE. Our generator injects a batch effect deliberately, but the
   pipeline does not remove it.
3. **Survival analysis.** Fit a Cox proportional-hazards model on the top
   biomarkers and validate the stratification on held-out cohorts. This is
   the point of a biomarker — predicting outcome, not just class.
4. **Nested cross-validation throughout.** Report performance from the outer
   CV of a properly nested design; the current `nested_cv_evaluate` is the
   starting point.
5. **Pre-registered thresholds.** Fix α, min |log₂FC|, and effect-size cutoffs
   *before* looking at data. Post-hoc threshold tuning is the most common
   source of irreproducible biomarker results.
6. **Wet-lab validation.** The ultimate test. No classifier result, however
   clean, replaces an assay.

None of these are implementable inside this repo. They are the reason a demo
like this is a starting point, not a study.

---

## Q11. Is this "pharmaceutical-grade" as the original description claimed?

**Short answer:** No, and the README no longer claims it.

**Long answer:** The original repo description used language like
"pharmaceutical-grade" and cited performance metrics (92.4% accuracy,
"450+ biomarkers") that came from a single synthetic run. Those claims have
been removed:

- The README's "Performance Metrics" section was replaced with an honest
  methodological description.
- The synthetic data was deliberately made harder (`effect_size = 0.5`) so
  models do not reach suspicious accuracy.
- The methods doc and this FAQ are written to be *auditable*, not persuasive.

"Pharmaceutical-grade" implies GxP-compliant processes, validated codebases,
audit trails, and regulatory review. A synthetic-data demo is not that. What
this repo *is*: a well-structured example of correct methodology for a
biomarker-discovery pipeline. That is a meaningfully different claim, and it
is the honest one.

---

## Q12. I see a `PROJECT_SUMMARY.md` file. What is it?

**Short answer:** A short changelog / summary of what changed from the
original version of the repo.

**Long answer:** `PROJECT_SUMMARY.md` documents the state of the codebase and
the historical fixes applied to it. It is intentionally short and factual.
For full methodological detail, use `docs/methods.md`; for figure descriptions,
use `docs/figures.md`; for the answers to skeptical questions, use this file.

The three docs split responsibilities:

| File | Audience | Content |
|---|---|---|
| `README.md` | anyone landing on the repo | what it is, how to run it |
| `docs/methods.md` | a reviewer auditing the science | formulas, thresholds, decisions |
| `docs/figures.md` | a reviewer reading the plots | what each figure shows, how to sanity-check |
| `docs/reviewer_faq.md` | a skeptical reviewer | this file |
| `PROJECT_SUMMARY.md` | the maintainer | what changed and why |

---

## Contact

If a question is not answered here and you would like it added, open an issue
or email f.e.amenaghawon@gmail.com. This FAQ is intended to be extended
whenever a new reviewer question appears.
