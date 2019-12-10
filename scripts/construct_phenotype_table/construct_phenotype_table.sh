# ARGS1: if you'd like to force running everything, put something here

# this is the script constructing phenotype table used for analysis
# Procedures:
# 1. Construct query
# 2. Query using ukbREST
# 3. Population assignment and QC
# 4. Phenotype QC


# construct query
echo '--> construct query: START'
if [[ ! -f ../../output/query_phenotypes.yaml ]] || [[ ! -z $1 ]]
then
  Rscript make_query.R \
    ../../misc/query_sample_filter_and_covar.yaml \
    ../../output/query_phenotypes.yaml
  echo '--> construct query: FINISHED'
else
  echo '--> construct query: output exists, SKIP'
fi
echo '##'

# query ukbREST
echo '--> query ukbREST: START'
if [[ ! -f ../../output/query_phenotypes_output.csv ]] || [[ ! -z $1 ]]
then
  curl -X POST \
    -H "Accept: text/csv" \
    -F file=@../../output/query_phenotypes.yaml \
    -F section=data \
    http://127.0.0.1:12345/ukbrest/api/v1.0/query > \
    ../../output/query_phenotypes_output.csv
  echo '--> query ukbREST: FINISHED'
else
  echo '--> query ukbREST: output exists, SKIP'
fi
echo '##'

# population assignment and QC
echo '--> population assignment and QC: START'
if [[ ! -f ../../output/query_phenotypes_with_population_qc.rds ]] || [[ ! -z $1 ]]
then
  Rscript pop_qc.R \
    ../../output/query_phenotypes_output.csv \
    ../../output/query_phenotypes_with_population_qc.rds
  echo '--> population assignment and QC: FINISHED'
else
  echo '--> population assignment and QC: output exists, SKIP'
fi
echo '##'

# phenotype QC
echo '--> phenotype QC: START'
if [[ ! -f ../../output/query_phenotypes_cleaned_up.csv ]] || [[ ! -z $1 ]]
then
  Rscript phenotype_qc.R \
    ../../output/query_phenotypes_with_population_qc.rds \
    ../../output/query_phenotypes_cleaned_up.csv
  echo '--> phenotype QC: FINISHED'
else 
  echo '--> phenotype QC: output exists, SKIP'
fi
echo '##'

# split British into training, validation, and test
echo '--> data split: START'
if [[ ! -d ../../output/data_split/ ]] || [[ ! -z $1 ]]
then
  Rscript data_split.R \
    ../../output/query_phenotypes_cleaned_up.csv \
    ../../output/data_split
  echo '--> data split: FINISHED'
else
  echo '--> data split: output exists, SKIP'
fi
echo '##'

