# Balanced MRD virtual-patient nested feature selection

 The virtual-patient table is `Data/PB_BM_merged_with_MRDflow_BCR.csv`. `Data/markov__deltaBIC2__final_subset_summary.csv` holds the original Markov winner scores used in the comparison figure. 

From MATLAB, change to this directory first.

```matlab
build_balanced_mrd_vp
collect_vp_mrd_mean_lambdas
run_balanced_mrd_nested_search
make_balanced_mrd_fold_overlay_figure
make_balanced_mrd_outer_auc_ba_figure
```

`build_balanced_mrd_vp` keeps 85 negative and 85 positive rows. Negative specimens contribute one virtual patient. Positive specimens contribute four, and nine of them contribute a fifth. Seed 1729.

`collect_vp_mrd_mean_lambdas` tunes ridge, lasso, and elastic net on all 10,400 eligible flow-MRD rows, with every row from a specimen kept in the same fold. The saved penalty is the mean of the 500 outer-training choices.

`run_balanced_mrd_nested_search` searches subsets with 10 outer and 10 inner folds and those fixed penalties. The fold overlay and the accuracy figure use each fold's Delta-BIC winner.
