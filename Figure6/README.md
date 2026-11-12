# Virtual-patient figures

This folder rebuilds the virtual-patient analyses from the files inside it. The fitting code is in `Functions`. The only cohort table is `Data/PB_BM_merged_with_MRDflow_BCR.csv`. Four small tables in `Data/` are the original-cohort results that the virtual-patient figures are compared with. The scripts write their own result tables into `outputs/` when they run.



Specimen-level out-of-fold performance, then the original-versus-virtual-patient figure:

```matlab
addpath('Benchmark')
compute_vp_mean_oof_bootstrap
compute_vp_cluster_bootstrap
make_vp_outer_oof_bootstrap_figure
```
