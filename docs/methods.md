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

