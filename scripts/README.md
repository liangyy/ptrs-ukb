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

* Run GWAS

```
$ screen -dmS gwas bash run-submit-and-loop_gwas_on_all_subsets.screen  # screen -dmS run_gwas bash run_gwas_on_all_subsets.screen
```
