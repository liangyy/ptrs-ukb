import argparse
parser = argparse.ArgumentParser(prog='gwas_on_subset.py', description='''
    Run GWAS (use dosage) given a list of phenotypes 
    which are preprocessed by `gwas_build_pheno_and_covar.py`
    and a list of individual subsets.
    Output the hail Table containing GWAS results. 
    The post-processing/formatting is done separately.
''')

## Genotype inputs
parser.add_argument('--bgen-path', required=True, help='''
    Template of bgen file, need to contain {chr_num} in place of chromosome number
''')
all_chrs = ','.join([ str(i) for i in range(1, 23) ])
parser.add_argument('--chrs', default=all_chrs, type=str, help='''
    The chromosome indexes separated by ',' to include. 
    If not specified, all autosomes will be used.
''')
parser.add_argument('--bgen-sample', required=True, help='''
    Sample file of bgen
''')
parser.add_argument('--bgen-index', default=None, help='''
    If not specified, assume the idx2 files are in the same directory
    and same name as bgen files.
    Otherwise, specify with {chr_num} as well
''')
parser.add_argument('--hail-block-size', type=int, default=128, help='''
    The block_size in hail.import_bgen
''')
parser.add_argument('--output-filename', required=True, help='''
    Filename of output (it is hail Table). 
    If the directory does not exist, it will be created.
    If extension is not .ht, .ht will be added.
''')

## Variant QC input
parser.add_argument('--variant-ht', help='''
    The hail Table containing variant QC of the input genotype.
''')
parser.add_argument('--variant-qc-yaml', help='''
    A YAML file specifying the filters for variant QC 
''')

## Phenotype/covariate input
parser.add_argument('--pheno-covar-path', help='''
    The path of hail Tables with {subset_name}
''')

## Subset inputs
parser.add_argument('--subset-yaml', help='''
    A YAML file specifying the name and path of individual ID lists.
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

# input sanity check
if '{chr_num}' not in args.bgen_path:
    logging.info('Wrong --bgen-path! It should contain {chr_num}! Exit')
    sys.exit()
if '{chr_num}' not in args.bgen_index:
    logging.info('Wrong --bgen-index! It should contain {chr_num}! Exit')
    sys.exit()
if '{subset_name}' not in args.pheno_covar_path:
    logging.info('Wrong --pheno-covar-path! It should contain {subset_name}! Exit')
    sys.exit()

# more on args
if args.bgen_index is None:
    args.bgen_index = args.bgen_path + '.idx2'

# some hail environment logging before run
logging.info('echo $PYSPARK_SUBMIT_ARGS')
os.system('echo $PYSPARK_SUBMIT_ARGS')

# initialize hail
logging.info('Initialize hail')
hl.init()

# load variant QC file and apply filter
logging.info('Start variant QC')
variant_qc_dic = helper.read_yaml(args.variant_qc_yaml)
logging.info('--> Read variant QC Table')
tstart = time.time()
variant_qc_all_ht = hl.read_table(args.variant_ht)
tend = time.time()
logging.info('--> Read variant QC Table FINISHED. {} seconds elapsed'.format(tend - tstart))
for i in list(variant_qc_dic.keys()):
    logging.info('--> Apply filter {}'.format(i))
    tstart = time.time()
    variant_qc_all_ht = helper.apply_variant_qc_filter(variant_qc_all_ht, i, variant_qc_dic[i])
    tend = time.time()
    logging.info('--> Apply filter {} FINISHED! {} seconds elapsed'.format(i, tend - tstart))
logging.info('Varint QC FINISHED.')

# load bgen
logging.info('Start to load genotype files with block_size = {}'.format(args.hail_block_size))
bgen_path = args.bgen_path.format(
    chr_num = '{' + args.chrs + '}'
)
bgen_sample = args.bgen_sample
bgen_idx_dict = {}
for i in range(1, 23):
    bgen = args.bgen_path.format(chr_num = i)
    bgen_idx = args.bgen_index.format(chr_num = i)
    bgen_idx_dict[bgen] = bgen_idx
tstart = time.time()
mt = hl.import_bgen(
    path = bgen_path, 
    sample_file = bgen_sample, 
    n_partitions = None, 
    index_file_map = bgen_idx_dict, 
    entry_fields = ['dosage'],  # I use dosage rather than GT
    block_size = args.hail_block_size,
    variants = variant_qc_all_ht
)
tend = time.time()
logging.info('Loading genotype FINISHED! {} seconds elapsed'.format(tend - tstart))

# count the number of variants being loaded
logging.info('Start to count the number of variants being loaded')
tstart = time.time()
nvariant = mt.count_rows()
tend = time.time()
logging.info('Counting the the number of variants FINISHED! nvariant = {} and {} seconds elapsed'.format(nvariant, tend - tstart))


# load the list of hail Tables as nested phenotypes and covariates 
logging.info('Start loading hail Table by individual subset')
subset_list_dic = helper.read_yaml(args.subset_yaml)
subset_ht_dic = {}
for subset_name in list(subset_list_dic.keys()):
    logging.info('--> Loading {}'.format(subset_name))
    tstart = time.time()
    subset_ht_dic[subset_name] = hl.read_table(args.pheno_covar_path.format(subset_name = subset_name))
    tend = time.time()
    logging.info('--> Loading {} FINISHED! {} seconds elapsed'.format(subset_name, tend - tstart))
subset_name = 'covariates'
logging.info('--> Loading {}'.format(subset_name))
tstart = time.time()
covar = hl.read_table(args.pheno_covar_path.format(subset_name = subset_name))
tend = time.time()
logging.info('--> Loading {} FINISHED! {} seconds elapsed'.format(subset_name, tend - tstart))
logging.info('Loading hail Table by individual subset FINISHED!')

# annotate genotype MatrixTable with these nested phenotypes and covariates
logging.info('Start annotating genotypes with phenotypes and covariates')
## do nested phenotype first
logging.info('--> Annotate with {}'.format('nested phenotypes'))
tstart = time.time()
annot_expr_ = {
    k : subset_ht_dic[k][mt.s] for k in list(subset_ht_dic.keys())
}
mt = mt.annotate_cols(**annot_expr_)
tend = time.time()
logging.info('--> Annotate with {} FINISHED! {} seconds elapsed'.format('nested phenotypes', tend - tstart))
## then do covariates
logging.info('--> Annotate with {}'.format('covariates'))
tstart = time.time()
mt = mt.annotate_cols(covariates = covar[mt.s])
tend = time.time()
logging.info('--> Annotate with {} FINISHED! {} seconds elapsed'.format('covariates', tend - tstart))

# prepare phenotypes and covariates into list of lists and list
logging.info('Start preparing `y` and `covariates` for `linear_regression_rows`')
pheno_list_of_lists = [ [ mt[i][j] for j in mt[i] ] for i in list(subset_ht_dic.keys()) ]
pheno_list_of_names = [ [ f'{i}_x_{j}' for j in mt[i] ] for i in list(subset_ht_dic.keys()) ]
covar_list = [ mt.covariates[i] for i in list(mt.covariates.keys()) ]
logging.info('Prepare `y` and `covariates` for `linear_regression_rows` FINISHED!')

# run GWAS
logging.info('Start running GWAS')
tstart = time.time()
gwas_out = hl.linear_regression_rows(
    y = pheno_list_of_lists,
    x = mt.dosage,
    covariates = [1] + covar_list,
    pass_through = ['varid', 'rsid']
)
gwas_out = gwas_out.annotate_globals(phenotypes = pheno_list_of_names)
tend = time.time()
logging.info('Running GWAS FINISHED! {} seconds elapsed'.format(tend - tstart))

# write GWAS results onto disk
logging.info('Start writing GWAS result to disk')
tstart = time.time()
## if target folder does not exist, create it
target_folder = os.path.dirname(args.output_filename)
if not os.path.exists(target_folder) and target_folder is not '':
    os.makedirs(target_folder)
## check if extension of output file is .ht, if not add it
filename, file_extension = os.path.splitext(args.output_filename)
if file_extension != 'ht':
    filename = filename + '.ht'
gwas_out.write(filename, overwrite = True)
tend = time.time()
logging.info('Writing GWAS result to disk FINISHED! {} seconds elapsed'.format(tend - tstart))
