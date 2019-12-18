## About

Here we perform GCTA analysis to estimate the variation of phenotype can be explained by predicted expression. 
We just frame it as "gcta regulability" which is not precise so it is just a name.

Here we will implement two analysis pipelines. 

1. Use GCTA software
2. Use hail

There is not strict preference, and we just want to make side-by-side comparison for fun.

To make life easier, we will use European subset 1 test individuals along with other population that passed our phenotype QC's.

## Input

1. Predicted expression stored in hdf5 format
2. Phenotype/covariate table and individual lists

Regarding the covariates, we will use the same set of covariates as GWAS does (20 PCs + 5 covariates on sex and age).

## Workflow & output

All pain is about formatting and things going with it.

1. Extract subset of predicted expression for the individuals of interest
2. Extract subset of phenotypes and covariates for the individuals of interest
3. Save 1 and 2 as TSV so that it is easily readable by hail
4. Prepare GRM in GZ format for gcta run and also prepare phenotype and covariate files as well
