import argparse
parser = argparse.ArgumentParser(prog='estimate_hsq_multi_tissue.py', description='''
    Fit mixed effect model Y ~ X + Z + e 
    using hail LinearMixedModel internally.
    Where allow a sequence of predictor matrice.
''')
parser.add_argument('--trait-table', required=True, help='''
    Trait table in TSV format.
    Specify the column name for individual ID
    by adding '::[colname]' after file path.
    If not specified, it assumes there is a 
    column called 'eid' as individual ID.
''')
parser.add_argument('--covar-table', required=True, help='''
    Covariate table in TSV format.
    Specify the column name for individual ID
    by adding '::[colname]' after file path.
    If not specified, it assumes there is a 
    column called 'eid' as individual ID.
    ALL other columns will be treated as 
    quantitative variables and will be included
    in the analysis.
''')
parser.add_argument('--predictor-table-list', required=True, help='''
    A list of predictor matrix in TSV format (filenames are separated by '::')
    (predictor in rows and individual in columns).
''')
parser.add_argument('--predictor-gene-column', required=True, help='''
    Specify the column name for predictor ID
    by adding '::[colname]' after file path.
    If not specified, it assumes no predictor ID.
''')
parser.add_argument('--output', required=True, help='''
    Output file name (TSV format)
''')
parser.add_argument('--inv-norm-predictor', action="store_true", help='''
    If you'd like to inverse normalize predictor
''')
parser.add_argument('--standardize-predictor', action="store_true", help='''
    If you'd like to standardize predictor
''')
parser.add_argument('--with-intercept', action="store_true", help='''
    If you'd like to add intercept as fixed effect
''')
parser.add_argument('--mode', default="naive", help='''
    mode for preparing gcta input: naive as default
''')
args = parser.parse_args()

import logging, time, sys
import numpy as np
import pandas as pd
import gcta_helper as ghelper
import hail as hl

def parse_input(tag, default):
    if '::' in tag:
        o = tag.split('::')
        return o[0], o[1]
    else:
        return tag, default

# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

# input sanity check
if args.inv_norm_predictor == True and args.standardize_predictor == True:
    logging.info('This script does not support doing inverse normalization and standardization simultaneously ..')
    logging.info('Exit!')
    sys.exit()

# read in trait table
logging.info('Loading trait table')
filename, trait_colname = parse_input(args.trait_table, 'eid')
df_trait = ghelper.tsv_to_pd_df(filename, indiv_col = trait_colname)
indiv_pool = np.array(df_trait[trait_colname].to_list())
logging.info('--> Current sample size = {}'.format(indiv_pool.shape[0]))

# read in covariate table
logging.info('Loading covariate table')
filename, covar_colname = parse_input(args.covar_table, 'eid')
df_covar = ghelper.tsv_to_pd_df(filename, indiv_col = covar_colname)
indiv_pool = np.intersect1d(indiv_pool, df_covar[covar_colname].to_list())
logging.info('--> Current sample size = {}'.format(indiv_pool.shape[0]))

# read in predictor matrix
logging.info('Loading predictor matrice. Looping over all inputs')
predictor_tables = args.predictor_table_list.split('::')
if args.mode == 'naive':
    piled_pred_expr = None
    ntotal = len(predictor_tables)
    for i in range(ntotal):
        filename = predictor_tables[i]
        colname = args.predictor_gene_column
        df_pred_expr = ghelper.tsv_to_pd_df(filename, indiv_col = colname)
        df_pred_expr = df_pred_expr.drop(columns = [colname])
        if piled_pred_expr is None:
            piled_pred_expr = df_pred_expr
        else:
            piled_pred_expr = pd.concat((piled_pred_expr, df_pred_expr), ignore_index = True)
        indiv_pool = np.intersect1d(indiv_pool, piled_pred_expr.columns.to_list())
        logging.info('--> MODE {}, {}/{}, Current sample size = {}'.format(args.mode, i + 1, ntotal, indiv_pool.shape[0]))
elif args.mode == 'tissue_svd':
    # loop over all genes and do svd for each gene.
    # first of all, collect pred expr and genes
    list_pred_expr = []
    col_order = None
    genes = set()
    for i in range(ntotal):
        filename = predictor_tables[i]
        colname = args.predictor_gene_column
        df_pred_expr = ghelper.tsv_to_pd_df(filename, indiv_col = colname)
        # df_pred_expr = df_pred_expr.drop(columns = [colname])
        if col_order is None:
            col_order = list(df_pred_expr.columns)
            gene_idx = col_order.index(colname)
            col_order.pop(gene_idx)
            col_order.append(colname)
        else:
            df_pred_expr = df_pred_expr[col_order]
        list_pred_expr.append(df_pred_expr)
        genes.union(set(df_pred_expr[colname].to_list()))
    # second, loop over gene and do svd
    piled_pred_expr = None
    for g in genes:
        _mat = []
        for ele in list_pred_expr:
            _tmp = list(ele.loc[ele[colname] == g])
            _tmp.pop(-1)
            _mat.append(_tmp)
        _mat = np.array(_mat)
        
        

# organize tables and matrix so that they match by individual ordering
logging.info('Organizing tables and matrix by individual ordering')
df_trait = df_trait[np.isin(df_trait[trait_colname], indiv_pool)].sort_values(by = trait_colname)
df_covar = df_covar[np.isin(df_covar[covar_colname], indiv_pool)].sort_values(by = covar_colname)
piled_pred_expr = piled_pred_expr[df_trait[trait_colname].to_list()]
df_trait = df_trait.drop(columns = [trait_colname])
df_covar = df_covar.drop(columns = [covar_colname])
trait_names = df_trait.columns

# convert to numpy darray
logging.info('Preparing input to model call')
np_trait = df_trait.to_numpy()
np_covar = df_covar.to_numpy()
np_pred_expr = piled_pred_expr.to_numpy()
logging.info('--> Removing constant predictors from matrix. # predictors = {}'.format(np_pred_expr.shape[0]))
np_pred_expr, _ = ghelper.remove_constant_row(np_pred_expr, np.ones(np_pred_expr.shape[0]))
logging.info('--> After removing constant predictors, # predictors = {}'.format(np_pred_expr.shape[0]))
if args.inv_norm_predictor == True:
    logging.info('--> Inverse normalizing predictors')
    np_pred_expr = ghelper.inv_norm_row(np_pred_expr)
if args.standardize_predictor is True:
    logging.info('--> Standardizing predictors')
    np_pred_expr = ghelper.standardize_row(np_pred_expr)
if args.with_intercept is True:
    logging.info('--> Adding intercept')
    np_covar = np.append(np.ones([np_covar.shape[0], 1]), np_covar, axis = 1)

# loop over traits and fit model
logging.info('Fitting models')
logging.info('--> Computing K matrix (Z Z^t)')
K_mat = np.dot(np_pred_expr.transpose(), np_pred_expr) / np_pred_expr.shape[0]
init_nan = np.empty(len(trait_names))
init_nan[:] = np.nan
out = pd.DataFrame({
    'trait': trait_names, 
    'num_predictors': np_pred_expr.shape[0], 
    'num_samples' : np_pred_expr.shape[1], 
    'h_sq': init_nan,
    'h_sq_se': init_nan
})
for i in range(len(trait_names)):
    logging.info('--> Working on {}'.format(trait_names[i]))
    grm, _ = hl.stats.LinearMixedModel.from_kinship(
        y = np_trait[:, i],
        x = np_covar,
        k = K_mat
    )
    try:
        grm.fit()
        out['h_sq'][i] = grm.h_sq
        out['h_sq_se'][i] = grm.h_sq_standard_error
    except:
        pass

# save results
logging.info('Saving results')
out.to_csv(args.output, index = None, sep = '\t')

