This directory holds the collection of scripts used for obtaining phenotypes, covariates, and population assignment for the T2D (using T2D and HbA1c in UKB) analysis.
To run the pipeline, do

```
$ bash construct_phenotype_table.sh 1  # the additional '1' means you force re-runing the pipeline, otherwise it will skip steps if the output has already existed
```
