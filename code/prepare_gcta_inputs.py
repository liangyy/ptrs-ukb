import argparse
parser = argparse.ArgumentParser(prog='prepare_gcta_inputs.py', description='''
    Convert predicted expression and phenotype/covariates table to GCTA inputs
''')

parser.add_argument('--pred-expr', required=True, help='''
    Predicted expression TSV
''')
parser.add_argument('--pheno', required=True, help='''
    Phenotype TSV
''')
parser.add_argument('--covar', required=True, help='''
    Covariate TSV
''')
parser.add_argument('--output-grm-prefix', required=True, help='''
    Output prefix of GRM
''')
parser.add_argument('--output-pheno-prefix', required=True, help='''
    Output prefix of phenotypes
''')
parser.add_argument('--output-covar', required=True, help='''
    Output filename of covariates
''')
parser.add_argument('--indiv-colname', default='eid', help='''
    Column name of individual ID in input
''')


args = parser.parse_args()

import pandas as pd
import logging, sys
import gcta_helper as ghelper

# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

logging.info('Read predicted expression')
pred_expr = pd.read_csv(args.pred_expr, sep = '\t', header = 0)
gene_in_pred_expr = pred_expr['gene'].to_list()
pred_expr_mat = pred_expr.drop(columns = ['gene'])

logging.info('Inverse normalize')
inv_norm_pred_expr_mat = ghelper.inv_norm_row(pred_expr_mat)

logging.info('Build GRM')
grm = ghelper.format_to_gcta_grm(inv_norm_pred_expr_mat)

logging.info('Output GRM GZ')
grm.to_csv(args.output_grm_prefix + '.grm.gz', header = None, index = None, sep = '\t', compression = 'gzip')
grm_indiv = pred_expr_mat.columns.to_list()
logging.info('Output GRM sample ID')
pd.DataFrame({
    1 : grm_indiv, 
    2 : grm_indiv
}).to_csv(
    args.output_grm_prefix + '.grm.id', 
    header = None, 
    index = None, 
    sep = '\t'
)

logging.info('Read phenotypes')
pheno = pd.read_csv(args.pheno, header = 0, sep = '\t')
colnames = pheno.columns
for c in colnames:
    if c != args.indiv_colname:
        pheno[[args.indiv_colname, args.indiv_colname, c]].to_csv(args.output_pheno_prefix + c + '.pheno', header = None, index = None, sep = '\t')

logging.info('Read covariates')
covar = pd.read_csv(args.covar, header = 0, sep = '\t')
colnames = covar.columns
colnames = [ c for c in colnames if c != args.indiv_colname ]  # remove all columns with name args.indiv_colname
covar[[args.indiv_colname, args.indiv_colname] + colnames].to_csv(args.output_covar, header = None, index = None, sep = '\t')

