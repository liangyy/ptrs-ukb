Here it generates a subset among all British individuals passing QC steps in `../construct_phenotype_table/`.
The output will be used to build the reference LD panel for LD-clumping step in PRS building.

To run:

```
$ Rscript get_subset.R \
  --pheno_table ../../output/query_phenotypes_cleaned_up.csv \
  --indiv_col eid \
  --pop_col meaning \
  --pop British \
  --nindiv 2000 \
  --output ../../output/indiv_list_as_ld_reference.txt
```
