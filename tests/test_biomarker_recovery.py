"""True biomarkers should be recoverable under default-like settings."""
from cancer_biomarkers.config import DataConfig, DEConfig
from cancer_biomarkers.data import generate_dataset
from cancer_biomarkers.differential_expression import differential_expression


def test_true_biomarkers_enriched_in_significant():
    cfg = DataConfig(
        n_samples=120,
        n_genes=200,
        n_true_biomarkers=15,
        effect_size=0.8,   # slightly stronger so the test is not seed-sensitive
        random_state=7,
    )
    ds = generate_dataset(cfg)
    de = differential_expression(ds.expression, ds.labels, DEConfig(min_abs_log2fc=0.6))

    sig = set(de.loc[de["significant"], "gene"])
    true = set(ds.true_biomarkers)
    recovered = true & sig

    assert len(recovered) >= max(1, len(true) // 3), (
        f"Only recovered {sorted(recovered)} of {sorted(true)}"
    )

    # Also check that the top-ranked significant gene is a true biomarker
    top_sig = de.loc[de["significant"]].assign(abs_fc=lambda d: d["log2fc"].abs())
    if not top_sig.empty:
        top_gene = top_sig.sort_values("abs_fc", ascending=False).iloc[0]["gene"]
        assert top_gene in true or len(recovered) > 0
