import argparse
parser = argparse.ArgumentParser(prog='subset_pheno_covar_by_indiv.py', description='''
    Input1: pheno_covar table obtained from ukbREST along with post-QCs
    Input2: YAML file defining which columns are phenotypes and covariates
    Input3: a list of individual ID
    Output: the phenotype and covariate for the subset of individuals
''')

parser.add_argument('--pheno-covar-csv', required=True, help='''
    Phenotype table obtained from ukbREST 
''')
parser.add_argument('--pheno-covar-yaml', required=True, help='''
    YAML file telling which columns are phenotype and covariates
''')
parser.add_argument('--indiv-list', required=True, help='''
    The list of individuals (it can have several columns but the first one 
    will be treated as individual ID)
''')
parser.add_argument('--output-pheno', required=True, help='''
    Phenotype table for subset individuals
''')
parser.add_argument('--output-covar', required=True, help='''
    Covariate table for subset individuals
''')
parser.add_argument('--indiv-colname', default='eid', help='''
    Column name of individual ID in input
''')

args = parser.parse_args()

import pandas as pd
import numpy as np
import h5py
import logging, os, time, sys
import my_hail_helper as hail_helper
import gwas_helper


# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

# load phenotypes and covariates (Exactly the same as gwas_build_pheno_and_covar.py)
logging.info('Start loading phenotypes and covariates (the full table)')
pheno_covar_dic = gwas_helper.read_yaml(args.pheno_covar_yaml)
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

# read individual list
logging.info('Read individual list')
indiv_list = hail_helper.read_indiv_list(args.indiv_list)

# subsetting
trait_sub = hail_helper.subset_by_col(trait, args.indiv_colname, indiv_list)
covar_sub = hail_helper.subset_by_col(covar, args.indiv_colname, indiv_list)

# save as TSV
trait_sub.to_csv(args.output_pheno, header = True, index = None, sep = '\t')
covar_sub.to_csv(args.output_covar, header = True, index = None, sep = '\t')
