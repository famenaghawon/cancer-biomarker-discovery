from cancer_biomarkers.config import DataConfig, MLConfig, PathConfig
from cancer_biomarkers.pipeline import run


def test_pipeline_runs_end_to_end(tmp_path):
    paths = PathConfig(results_dir=tmp_path, figures_dir=tmp_path / "figures")

    summary = run(
        data_cfg=DataConfig(n_samples=40, n_genes=60, n_true_biomarkers=5),
        ml_cfg=MLConfig(test_size=0.25, n_splits=3, top_k_features=10),
        paths=paths,
    )

    # Output files
    assert (tmp_path / "analysis_summary.json").exists()
    assert (tmp_path / "expression_full.parquet").exists()
    assert (tmp_path / "significant_biomarkers.csv").exists()
    assert (tmp_path / "figures" / "volcano.png").exists()
    assert (tmp_path / "figures" / "heatmap.png").exists()
    assert (tmp_path / "figures" / "roc_curve.png").exists()

    # Summary schema
    assert "models" in summary["ml"]
    assert set(summary["ml"]["models"]) == {"logistic_regression", "random_forest"}
    assert "top_k_recovery" in summary
    recovery = summary["top_k_recovery"]
    assert recovery["k"] > 0
    assert 0.0 <= recovery["recovery_rate"] <= 1.0

    # Nested CV should not be present when the env var is unset
    assert "ml_nested" not in summary

    # Metrics in range
    for name, m in summary["ml"]["models"].items():
        for k in ("accuracy", "precision", "recall", "f1"):
            assert 0.0 <= m[k] <= 1.0, f"{name}.{k} out of range"
        assert 0.0 <= m["cv_roc_auc_mean"] <= 1.0
