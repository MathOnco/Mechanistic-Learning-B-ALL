# Nested feature selection

This folder selects Markov and flow subsets, chooses the Delta-BIC winner in each outer fold, and scores that winner on the held-out fold. The only input tables are `Data/markov.csv` and `Data/flow.csv`. The scripts write their own result tables when they run.

From MATLAB, change to this directory first.

```matlab
run_exhaustive_nested_10fold
run_flow_abundance_nested_10fold
select_final_subset_nested_delta_bic(fullfile(pwd,'outputs_10fold'), "markov")
select_final_subset_nested_delta_bic(fullfile(pwd,'outputs_10fold'), "flow")
addpath('Figures')
make_deltaBIC_fold_overlay_figure
make_deltaBIC_selected_auc_ba_comparison_figure
```

The first two calls are the 10-by-10 nested search. The penalty is fixed at the mean lambda from the classifier benchmark. `select_final_subset_nested_delta_bic` picks each fold's Delta-BIC winner and scores that winner on the outer-test rows. The two figure functions read those tables and write `outputs/figures/`.

A short installation check, which does not run the full search:

```matlab
setenv('EXHAUSTIVE_SMOKE_MODE','1')
run_exhaustive_nested_10fold
```
