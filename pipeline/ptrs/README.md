Given predicted expression matrix and S-PrediXcan gene-level effect size.
Compute PTRS for a subset of individuals.

The procedure is:

1. Given S-PrediXcan results, define PTRS models
    - Exclude genes with huge effect size
    - Construct PTRS model on the basis of p-value cutoff
2. Apply PTRS model to the subset individuals

