# Classification Tasks

The pipeline supports three classification tasks with clearly different difficulty levels. Each one teaches something different about what is and is not predictable from expression data alone.

Pass the `--task` flag to `scripts/run_pipeline.py` to choose:

```bash
python scripts/run_pipeline.py --data tcga_brca --task tumor_vs_normal
python scripts/run_pipeline.py --data tcga_brca --task subtype
python scripts/run_pipeline.py --data tcga_brca --task stage
