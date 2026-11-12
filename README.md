# Mechanistic learning figure code

MATLAB scripts and input tables used to reproduce the classifier and virtual-patient figures from the mechanistic-learning paper.

## Requirements

- MATLAB (R2021a or later recommended)
- Statistics and Machine Learning Toolbox

Each figure folder is self-contained. Change into that folder in MATLAB before running the commands in its README.

## Contents

| Folder | What it reproduces |
|---|---|
| `Figure4/` | Nested classifier benchmark (ridge, lasso, elastic net, random forest) |
| `Figure5/` | Nested exhaustive feature selection and Delta-BIC figures |
| `Figure6/` | Virtual-patient versus original-cohort performance figures |
| `Figure7/` | Balanced flow-MRD virtual-patient nested feature selection |
| `Supplementary data/` | Summary CSV tables for the published figures |
| `Supplementary figures/` | Scripts and PNG panels for the supplementary figures |


## Notes

- Paths are resolved from each script location; no machine-specific directories are required.
- Generated figures are written as PNG files under each folder's `outputs/` directory.
