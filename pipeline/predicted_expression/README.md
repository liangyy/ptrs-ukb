Here we build pipeline starting from 

1. UKB bgen files 
2. A PrediXcan db
3. A list of individuals

and compute the predicted expression on that PrediXcan db.

The workflow is 

1. Subset the original bgen file on both individuals in the list and variants in PrediXcan db
2. Format the subsetted genotype file
3. Run `PrediXcan.py`
