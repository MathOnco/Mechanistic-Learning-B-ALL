# Classifier benchmark

This folder reproduces the nested cross-validation comparison of ridge, lasso, elastic net, and random forest. Predictor views are flow proportion, Markov transitions, clinical/genomic variables, and Markov plus clinical/genomic. Outcomes are BCR::ABL1 and flow MRD.

Everything required is in this folder. From MATLAB, change to this directory first.

Redraw the figure from the saved intervals:

```matlab
recreate_classifier_benchmark_figure
```

That reads `benchmark_all_bootstrap_95CI.csv` and writes `outputs/Figure_R1_classifier_benchmark.png` and `.pdf`.

Refit the benchmark:

```matlab
analysis_01_nested_classifier_benchmark
```

The refit uses `Data/markov.csv` and `Data/flow.csv`, 100 repeats of 5-by-5 nested cross-validation, and 2,000 bootstrap resamples. New tables are written to `outputs/`. A short installation check is:

```matlab
setenv('BENCHMARK_SMOKE_MODE','1')
analysis_01_nested_classifier_benchmark
```
