# ARGS1: if you'd like to force running everything, put something here

# this is the script constructing phenotype table used for analysis
# Procedures:
# 1. New population assignment and QC
# 2. Phenotype QC


# population assignment and QC
echo '--> population assignment and QC: START'
if [[ ! -f ../../output/new2_query_phenotypes_with_population_qc.rds ]] || [[ ! -z $1 ]]
then
  Rscript new_pop_qc.R \
    ../../output/query_phenotypes_output.csv \
    ../../output/new2_query_phenotypes_with_population_qc.rds
  echo '--> population assignment and QC: FINISHED'
else
  echo '--> population assignment and QC: output exists, SKIP'
fi
echo '##'

# phenotype QC
echo '--> phenotype QC: START'
if [[ ! -f ../../output/new2_query_phenotypes_cleaned_up.csv ]] || [[ ! -z $1 ]]
then
  Rscript ../construct_phenotype_table/phenotype_qc.R \
    ../../output/new2_query_phenotypes_with_population_qc.rds \
    ../../output/new2_query_phenotypes_cleaned_up.csv
  echo '--> phenotype QC: FINISHED'
else 
  echo '--> phenotype QC: output exists, SKIP'
fi
echo '##'

# split British into training, validation, and test
echo '--> data split: START'
if [[ ! -d ../../output/new2_data_split/ ]] || [[ ! -z $1 ]]
then
  Rscript new_data_split.R \
    ../../output/new2_query_phenotypes_cleaned_up.csv \
    ../../output/new2_data_split
  echo '--> data split: FINISHED'
else
  echo '--> data split: output exists, SKIP'
fi
echo '##'

