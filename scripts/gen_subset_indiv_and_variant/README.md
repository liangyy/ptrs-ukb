Here it generates a subset among all British individuals passing QC steps in `../construct_phenotype_table/`.
The output will be used to build the reference LD panel for LD-clumping step in PRS building.

To get subset individuals:

```
$ Rscript get_subset_indiv.R \
  --pheno_table ../../output/query_phenotypes_cleaned_up.csv \
  --indiv_col eid \
  --pop_col meaning \
  --pop British \
  --nindiv 2000 \
  --output ../../output/indiv_list_as_ld_reference.txt
```

To get the list of variants

```
$ cat /vol/bmd/yanyul/UKB/gwas_on_subset/gwas_runs_in_tsv/gwas_in_tsv_subset13_x_height.tsv| cut -f 12 | tail -n +2 > ../../output/variant_list_as_ld_reference.txt 
```

