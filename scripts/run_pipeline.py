"""CLI entry point: `python scripts/run_pipeline.py`."""
from __future__ import annotations

import sys
from pathlib import Path

# Make `src/` importable when run as a script from repo root
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cancer_biomarkers.pipeline import run  # noqa: E402

if __name__ == "__main__":
    run()
