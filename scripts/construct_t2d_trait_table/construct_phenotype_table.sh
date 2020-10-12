# ARGS1: if you'd like to force running everything, put something here

# this is the script constructing phenotype table used for analysis
# Procedures:
# 1. Construct query
# 2. Query using ukbREST
# 3. Population assignment and QC
# 4. Phenotype QC


# construct query
echo '--> construct query: START'
if [[ ! -f ../../output/query_t2d.yaml ]] || [[ ! -z $1 ]]
then
  Rscript make_query.R \
    ../../misc/query_sample_filter_and_covar.yaml \
    ../../output/query_t2d.yaml
  echo '--> construct query: FINISHED'
else
  echo '--> construct query: output exists, SKIP'
fi
echo '##'

# DONE by Sabrina Mi
# # query ukbREST
# echo '--> query ukbREST: START'
# if [[ ! -f ../../output/t2d_query.csv ]] || [[ ! -z $1 ]]
# then
#   curl -X POST \
#     -H "Accept: text/csv" \
#     -F file=@../../output/query_t2d.yaml \
#     -F section=data \
#     http://127.0.0.1:12345/ukbrest/api/v1.0/query > \
#     ../../output/t2d_query.csv
#   echo '--> query ukbREST: FINISHED'
# else
#   echo '--> query ukbREST: output exists, SKIP'
# fi
# echo '##'

# phenotype QC
echo '--> phenotype QC: START'
if [[ ! -f ../../output/t2d_query_cleaned_up.csv ]] || [[ ! -z $1 ]]
then
  Rscript phenotype_qc.R \
    ../../output/t2d_query.csv \
    ../../output/t2d_query_cleaned_up.csv
  echo '--> phenotype QC: FINISHED'
else 
  echo '--> phenotype QC: output exists, SKIP'
fi
echo '##'

