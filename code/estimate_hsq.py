import argparse
parser = argparse.ArgumentParser(prog='estimate_hsq.py', description='''
    Fit mixed effect model Y ~ X + Z + e 
    using hail LinearMixedModel internally.
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
parser.add_argument('--predictor-table', default=None, help='''
    Predictor matrix in TSV format 
    (predictor in rows and individual in columns).
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
parser.add_argument('--gene-list', default=None, help='''
    If you'd like to limit the REML analysis to a specific set of genes, 
    specify it here.
''')
parser.add_argument('--grm-gz', default=None, help='''
    If you'd like to use GRM matrix, set it here. --predictor-table will be
    high priority if both are specified (will report min, 25%, 50%, 75%, max N predictors).
''')
args = parser.parse_args()

import logging, time, sys
import gzip
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

def load_list(fn):
    o = []
    with open(fn, 'r') as f:
        for i in f:
            i = i.strip()
            o.append(i)
    return o

def load_grm_gz(prefix):
    fgrm = prefix + '.grm.gz'
    igrm = prefix + '.grm.id'
    ids = pd.read_csv(igrm, sep = '\t', header = None)
    indivs = list(ids.iloc[:, 0].astype(str).values)
    grm_mat = np.zeros((len(indivs), len(indivs)))
    grm_mat = grm_mat * np.float('nan')
    n_list = []
    with gzip.open(fgrm, 'rt') as f:
        for i in f:
            i = i.strip().split('\t')
            ii = int(i[0]) - 1
            jj = int(i[1]) - 1
            nn = int(float(i[2]))
            vv = float(i[3])
            grm_mat[ii, jj] = vv
            grm_mat[jj, ii] = vv
            n_list.append(nn)
    return grm_mat, indivs, np.array(n_list)

def _summarize_grm_N(nn):
    qq = np.quantile(nn, [0, 0.25, 0.5, 0.75, 1])
    fq = [ "{:.1f}".format(i) for i in list(qq) ]
    return ','.join(fq)

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

# read gene list if there is gene list
gene_list = None
if args.gene_list is not None:
    gene_list = load_list(args.gene_list)

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

if args.predictor_table is not None:
    # read in predictor matrix
    logging.info('Loading predictor matrix')
    filename, colname = parse_input(args.predictor_table, '')
    df_pred_expr = ghelper.tsv_to_pd_df(filename, indiv_col = colname)
    if gene_list is not None:
        df_pred_expr = df_pred_expr[ df_pred_expr.gene.isin(gene_list) ].reset_index(drop=True)
    df_pred_expr = df_pred_expr.drop(columns = [colname])
    indiv_pool = np.intersect1d(indiv_pool, df_pred_expr.columns.to_list())
    logging.info('--> Current sample size = {}'.format(indiv_pool.shape[0]))
    logging.info('--> Current number of genes = {}'.format(df_pred_expr.shape[0] - 1))
elif args.grm_gz is not None:
    # load GRM GZ
    logging.info('Loading GRM GZ matrix')
    K_mat, indiv_list, N_grm = load_grm_gz(args.grm_gz)
    indiv_pool = np.intersect1d(indiv_pool, indiv_list)
    logging.info('--> Current sample size = {}'.format(indiv_pool.shape[0]))
else:
    raise ValueError('No predictor being input.')

# organize tables and matrix so that they match by individual ordering
logging.info('Organizing tables and matrix by individual ordering')
df_trait = df_trait[np.isin(df_trait[trait_colname], indiv_pool)].sort_values(by = trait_colname)
df_covar = df_covar[np.isin(df_covar[covar_colname], indiv_pool)].sort_values(by = covar_colname)
if args.predictor_table is not None:
    df_pred_expr = df_pred_expr[df_trait[trait_colname].to_list()]
elif args.grm_gz is not None:
    df_idx_grm = pd.DataFrame({'idx': np.arange(K_mat.shape[0]), 'indiv': indiv_list})
    df_idx_target = pd.DataFrame({'indiv': df_trait[trait_colname].astype(str)})
    df_idx_target = pd.merge(df_idx_target, df_idx_grm, on = 'indiv')
    idx = df_idx_target.idx.values
    K_mat = K_mat[idx, :][:, idx]
df_trait = df_trait.drop(columns = [trait_colname])
df_covar = df_covar.drop(columns = [covar_colname])
trait_names = df_trait.columns

# convert to numpy darray
logging.info('Preparing input to model call')
np_trait = df_trait.to_numpy()
np_covar = df_covar.to_numpy()
if args.predictor_table is not None:
    np_pred_expr = df_pred_expr.to_numpy()
    logging.info('--> Removing constant predictors from matrix. # predictors = {}'.format(np_pred_expr.shape[0]))
    np_pred_expr, _ = ghelper.remove_constant_row(np_pred_expr, np.ones(np_pred_expr.shape[0]))
    logging.info('--> After removing constant predictors, # predictors = {}'.format(np_pred_expr.shape[0]))
    if args.inv_norm_predictor == True:
        logging.info('--> Inverse normalizing predictors')
        np_pred_expr = ghelper.inv_norm_row(np_pred_expr)
    if args.standardize_predictor is True:
        logging.info('--> Standardizing predictors')
        np_pred_expr = ghelper.standardize_row(np_pred_expr)   
    logging.info('--> Computing K matrix (Z Z^t)')
    K_mat = np.dot(np_pred_expr.transpose(), np_pred_expr) / np_pred_expr.shape[0] 
    num_predictors = np_pred_expr.shape[0]
    num_samples = np_pred_expr.shape[1]
elif args.grm_gz is not None:
    num_samples = K_mat.shape[0]
    num_predictors = _summarize_grm_N(N_grm)
if args.with_intercept is True:
    logging.info('--> Adding intercept')
    np_covar = np.append(np.ones([np_covar.shape[0], 1]), np_covar, axis = 1)
    
# loop over traits and fit model
logging.info('Fitting models')
init_nan = np.empty(len(trait_names))
init_nan[:] = np.nan
out = pd.DataFrame({
    'trait': trait_names, 
    'num_predictors': num_predictors, 
    'num_samples' : num_samples, 
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

