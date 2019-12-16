import argparse
parser = argparse.ArgumentParser(prog='gwas_build_pheno_and_covar.py', description='''
    Prepare the nested phenotypes and covariates for GWAS run.
''')


## Phenotype/covariate input
parser.add_argument('--pheno-covar-csv', help='''
    The CSV file queried from ukbREST 
    (there could be some extra filtering and QC as you see fit)
''')
parser.add_argument('--pheno-covar-yaml', help='''
    A YAML file specifying the path to full table and 
    the column names of phenotypes and covariates.
''')

## Subset inputs
parser.add_argument('--subset-yaml', help='''
    A YAML file specifying the name and path of individual ID lists.
''')

## Output prefix
parser.add_argument('--output-prefix', required=True, help='''
    Prefix of output file. 
    It will be [output-prefix].[subset_name].ht for subset phenotypes.
    It will be [output-prefix].covariates.ht for covariates
''')


args = parser.parse_args()

import hail as hl
import pandas as pd
import numpy as np
import logging, os, time, sys
import my_hail_helper as hail_helper
import gwas_helper as helper


# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

# some hail environment logging before run
logging.info('echo $PYSPARK_SUBMIT_ARGS')
os.system('echo $PYSPARK_SUBMIT_ARGS')

# initialize hail
logging.info('Initialize hail')
hl.init()

# load phenotypes and covariates
logging.info('Start loading phenotypes and covariates (the full table)')
pheno_covar_dic = helper.read_yaml(args.pheno_covar_yaml)
covar_names = pheno_covar_dic['covar_names']  # 'age_recruitment,sex,pc1,pc2'
pheno_names = pheno_covar_dic['pheno_names']  # 'ht,mcv,mch'
indiv_id = pheno_covar_dic['indiv_id']  # 'eid'
int_names = pheno_covar_dic['int_names']  # 'age_recruitment,sex'
str_names = pheno_covar_dic['str_names']  # 'eid'
logging.info('--> Read in CSV file as data.frame')
tstart = time.time()
covar, trait = hail_helper.read_and_split_phenotype_csv(
    args.pheno_covar_csv,
    pheno_names = pheno_names,
    covar_names = covar_names,
    indiv_id = indiv_id,
    int_names = int_names,
    str_names = str_names
)
tend = time.time()
logging.info('--> Read in CSV file as data.frame FINISHED! {} seconds elapsed'.format(tend - tstart))
# individual ID (the column key) is named 's' in genotype hail MatrixTable. 
# So that we change the column for indiv_id in pheno_covar table to 's' to make it consistent 
# (make life easier, sigh ...)
covar = covar.rename(columns = {indiv_id: 's'})
trait = trait.rename(columns = {indiv_id: 's'})
logging.info('--> Convert covariate data.frame to hail Table')
tstart = time.time()
covar = hail_helper.df_to_ht(covar, 's')
tend = time.time()
logging.info('--> Convert covariate data.frame to hail Table FINISHED! {} seconds elapsed'.format(tend - tstart))
logging.info('Loading phenotypes and covariates (the full table) FINISHED!')

# load the list of subsets and construct hail Tables 
logging.info('Start constructing subsets by individual')
subset_list_dic = helper.read_yaml(args.subset_yaml)
subset_ht_dic = {}
for subset_name in list(subset_list_dic.keys()):
    logging.info('--> Construct {}'.format(subset_name))
    tstart = time.time()
    subset_indiv_list = hail_helper.read_indiv_list(subset_list_dic[subset_name])
    sub_trait = hail_helper.subset_by_col(trait, 's', subset_indiv_list)
    sub_trait = hail_helper.df_to_ht(sub_trait, 's')  
    subset_ht_dic[subset_name] = sub_trait
    tend = time.time()
    logging.info('--> Construct {} FINISHED! {} seconds elapsed'.format(subset_name, tend - tstart))
logging.info('Constructing subsets by individual FINISHED!')
 
# write phenotypes and covariates onto disk
logging.info('Start writing phenotype and covariates hail Table to disk')
## if target folder does not exist, create it
target_folder = os.path.dirname(args.output_prefix)
if not os.path.exists(target_folder) and target_folder is not '':
    os.makedirs(target_folder)
tstart = time.time()
for subset_name in list(subset_ht_dic.keys()):
    logging.info('--> Writing {} to disk'.format(subset_name))
    tstart = time.time()
    subset_ht_dic[subset_name].write('{out_prefix}.{subset_name}.ht'.format(out_prefix = args.output_prefix, subset_name = subset_name), overwrite = True)
    tend = time.time()
    logging.info('--> Writing {} to disk FINISHED! {} seconds elapsed'.format(subset_name, tend - tstart))
subset_name = 'covariates'
logging.info('--> Writing {} to disk'.format(subset_name))
tstart = time.time()
covar.write('{out_prefix}.{subset_name}.ht'.format(out_prefix = args.output_prefix, subset_name = subset_name), overwrite = True)
tend = time.time()
logging.info('--> Writing {} to disk FINISHED! {} seconds elapsed'.format(subset_name, tend - tstart))
logging.info('Writing phenotype and covariates hail Table to disk FINISHED! {} seconds elapsed'.format(tend - tstart))
