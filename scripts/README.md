Here are the collection of scripts for this and that, mostly about one-shot formatting and etc.

* Indexing UKB genotypes in BGEN format for HAIL: `index_ukb_bgen_for_hail.sh`

```
$ nohup bash index_ukb_bgen_for_hail.sh 2> index_ukb_bgen_for_hail.log &
```

* Transfer data from CRI to nuleus

Here, because of ssh connection issue on CRI, we take a detour: CRI -> local -> nucleus.
It relies on two scripts: 

1. `transfer_pred_models_from_cri_to_local.sh`
2. `transfer_pred_models_from_local_to_nucleus.sh`

On local machine, run:

```
$ screen -dmS cri_to_local bash transfer_pred_models_from_cri_to_local.sh
```

After it is done, run:

```
$ screen -dmS local_to_nucleus bash transfer_pred_models_from_local_to_nucleus.sh
```

* Construct phenotype table used in the analysis

See `construct_phenotype_table/`.

* Do variant QC

```
$ screen -dmS variant_qc bash run_variant_qc.screen
```

* Synchronize `../output` at local and nucleus

```
$ bash sync_output_on_nucleus.sh
``` 

* Test run of `../code/gwas_on_subset.py` 

```
$ screen -dmS test_gwas bash test_run_gwas_on_subset.screen
```

* Run GWAS (DEPRECATED)

This is deprecated.
The actual runs were done on GCP.
See submission scripts at `../submission_scripts/for_gcp/`

```
$ screen -dmS gwas bash run-submit-and-loop_gwas_on_all_subsets.screen  # screen -dmS run_gwas bash run_gwas_on_all_subsets.screen
```

* Run PRS

Pre-subset genotypes.

```
$ screen -dmS presubset_1 bash run_prs_pre_subset.screen /vol/bmd/yanyul/GitHub/ptrs-ukb/misc/prs_full_1st_quar.yaml presubset_1 
$ screen -dmS presubset_2 bash run_prs_pre_subset.screen /vol/bmd/yanyul/GitHub/ptrs-ukb/misc/prs_full_2nd_quar.yaml presubset_2
$ screen -dmS presubset_3 bash run_prs_pre_subset.screen /vol/bmd/yanyul/GitHub/ptrs-ukb/misc/prs_full_3rd_quar.yaml presubset_3
$ screen -dmS presubset_4 bash run_prs_pre_subset.screen /vol/bmd/yanyul/GitHub/ptrs-ukb/misc/prs_full_4th_quar.yaml presubset_4
```

Calculate PRS

```
$ screen -dmS skip_1 bash run_prs_skip_subset.screen /vol/bmd/yanyul/GitHub/ptrs-ukb/misc/prs_full_1st_quar.yaml presubset_1
```

