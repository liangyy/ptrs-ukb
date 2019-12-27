# About

The script to calculate PRS is at `../../code/ld_clump_prs.py`.
It takes a set of GWAS results along with LD-clumping and the target individual lists

1. Subset genotypes data by variants and individuals
2. Loop over GWAS results and compute PRS 

Step 1 and 2 could be done separately by specifying `--mode` in the script ARGS: 

1. `full` runs all altogether
2. `pre_subset` runs subsetting only
3. `skip_subset` assumes subsetting has been done and runs PRS calculation only

Step 1 has been done on `nucleus`.
See scripts at `../../scripts/run_prs_pre_subset.screen`.
I tried PRS calculation (step 2) on `nucleus` as well at `../../scripts/run_prs_skip_subset.screen` but for some reason, it did not end up well. 
The issues are:

* Even though the genotype has been pre-subsetted, the genome-wide calculation is still intensive on `nucleus`
* The script needs high I/O which is a bit heavy on `nucleus` 
* Disk space could be a factor since I was keeping many computationally intensive intermediate files

# Run job on GCP

Test run to estimate resources and running time.

```
```



