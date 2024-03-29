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
parser.add_argument('--predictor-table-prefix', required=True, help='''
    To use --predictor-table-yaml, this should be specified.
''')
parser.add_argument('--predictor-table-suffix', required=True, help='''
    To use --predictor-table-yaml, this should be specified.
''')
parser.add_argument('--predictor-table-yaml', required=True, help='''
    A YAML file containing the list of predictor matrix in TSV format 
    along with the name
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
    If you'd like to inverse normalize predictor.
    It will get ignored when --mode=tissue_svd_train
''')
parser.add_argument('--standardize-predictor', action="store_true", help='''
    If you'd like to standardize predictor.
    It will get ignored when --mode=tissue_svd_train
''')
parser.add_argument('--with-intercept', action="store_true", help='''
    If you'd like to add intercept as fixed effect
''')
parser.add_argument('--mode', default="naive", help='''
    mode for preparing gcta input: naive as default
    For --mode=tissue_svd_train, it will loop over all genes and perform evd on cor(T).
    And keep the eigen vectors with lambda / lambda_max > 1 / 30 and output would be 
    the PC's for each gene.
    For --mode=tissue_svd, it takes the output generated by --mode=tissue_svd_train.
    And calculate the projected value of T (standardized by (obs - mean) / std) 
    as T %*% [PC1, ..., PCk].
    And perform reml on the projected matrix
''')
parser.add_argument('--pc-model', default=None, help='''
    If --mode=tissue_svd, specify the path to the PCA results (generated by --mode=tissue_svd_train)
    Ideally, the --predictor-table-yaml should be kept the same for --mode=tissue_svd 
    with PC models from --mode=tissue_svd_train to ensure the same tissue set and the same ordering.
''')
args = parser.parse_args()

import logging, time, sys
import numpy as np
import pandas as pd
import gcta_helper as ghelper
import hail as hl
import gwas_helper
import h5py

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
if args.mode == 'tissue_svd' and args.pc_model is None:
    logging.info('If --mode={}, then --pc-model should be specified'.format(args.mode))
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
predictor_tables_dic = gwas_helper.read_yaml(args.predictor_table_yaml)  # args.predictor_table_list.split('::')
print(predictor_tables_dic)
predictor_table_keys = list(predictor_tables_dic.keys())
ntotal = len(predictor_table_keys)
if args.mode == 'naive':
    piled_pred_expr = None
    ntotal = len(predictor_table_keys)
    for i in range(ntotal):
        filename = args.predictor_table_prefix + predictor_tables_dic[predictor_table_keys[i]] + args.predictor_table_suffix
        colname = args.predictor_gene_column
        df_pred_expr = ghelper.tsv_to_pd_df(filename, indiv_col = colname)
        df_pred_expr = df_pred_expr.drop(columns = [colname])
        if piled_pred_expr is None:
            piled_pred_expr = df_pred_expr
        else:
            piled_pred_expr = pd.concat((piled_pred_expr, df_pred_expr), ignore_index = True)
        indiv_pool = np.intersect1d(indiv_pool, piled_pred_expr.columns.to_list())
        logging.info('--> MODE {}, {}/{}, Current sample size = {}'.format(args.mode, i + 1, ntotal, indiv_pool.shape[0]))
elif args.mode == 'tissue_svd_train' or args.mode == 'tissue_svd':
    # loop over all genes and do svd for each gene.
    # first of all, collect pred expr and genes
    list_pred_expr = {}
    col_order = None
    genes = set()
    for i in range(ntotal):
        logging.info('--> MODE {}, read in tissue {}/{}'.format(args.mode, i + 1, ntotal))
        filename = args.predictor_table_prefix + predictor_tables_dic[predictor_table_keys[i]] + args.predictor_table_suffix
        colname = args.predictor_gene_column
        df_pred_expr = ghelper.tsv_to_pd_df(filename, indiv_col = colname)
        # df_pred_expr = df_pred_expr.drop(columns = [colname])
        if col_order is None:
            col_order = list(df_pred_expr.columns)
            gene_idx = col_order.index(colname)
            col_order.pop(gene_idx)
            # col_order_no_gene = col_order.copy()
            col_order.append(colname)
        else:
            df_pred_expr = df_pred_expr[col_order]
        list_pred_expr[predictor_table_keys[i]] = df_pred_expr
        genes = genes.union(set(df_pred_expr[colname].to_list()))
    # second, loop over gene and do svd
    if args.mode == 'tissue_svd':
        piled_pred_expr = []
        with h5py.File(args.pc_model, 'r') as model_handle:
            avail_genes = set(model_handle.keys())
            counter = 0
            ngene = len(genes)
            check_n = int(ngene / 50)
            for g in genes:
                counter += 1
                if counter % check_n == 1 or ngene == counter:
                    logging.info('----> MODE {}, working on gene {}/{}'.format(args.mode, counter, ngene))
                if g not in avail_genes:
                    continue
                g_tissues = model_handle['{}/tissues'.format(g)][...].astype(str)
                _mat = []
                for ele_key in g_tissues:  # list_pred_expr.keys():
                    if ele_key not in list_pred_expr:
                        raise ValueError('key {} is not in input pred expr tables'.format(ele_key))
                    ele = list_pred_expr[ele_key]
                    _tmp = ele.loc[ele[colname] == g].to_numpy()
                    if _tmp.size != 0:
                        _tmp = _tmp[0, :][: (_tmp.shape[1] - 1)]
                        _mat.append(_tmp)
                    else:
                        continue
                _mat = np.array(_mat)
                if _mat.size != 0:
                    _mat = ghelper.standardize_row(_mat)
                    v =  model_handle['{}/v'.format(g)][...]
                    _p_mat = np.matmul(v.T, _mat)
                    piled_pred_expr.append(_p_mat)
        col_order.pop(-1)
        piled_pred_expr = pd.DataFrame(np.concatenate(piled_pred_expr).astype(float), columns = col_order)
    elif args.mode == 'tissue_svd_train':
        ngene = len(genes)
        check_n = int(ngene / 50)
        counter = 0
        with h5py.File(args.output, 'w') as out_handle:
            out_handle.create_dataset('tissue_list_in_order', data = np.array(predictor_table_keys).astype('S'))
            for g in genes:
                counter += 1
                if counter % check_n == 1 or ngene == counter:
                    logging.info('----> MODE {}, working on gene {}/{}'.format(args.mode, counter, ngene))
                _mat = []
                _tissues = []
                for ele_key in list_pred_expr.keys():
                    ele = list_pred_expr[ele_key]
                    _tmp = ele.loc[ele[colname] == g].to_numpy()
                    if _tmp.size != 0:
                        _tmp = _tmp[0, :][: (_tmp.shape[1] - 1)]
                        _mat.append(_tmp)
                        _tissues.append(ele_key)
                    else:
                        continue
                _mat = np.array(_mat)
                if _mat.size != 0:
                    _mat = ghelper.standardize_row(_mat)
                    TtT = np.matmul(_mat, _mat.T) / _mat.shape[1]
                    w, v = np.linalg.eigh(TtT.astype(float))
                    w, v = ghelper.truncate_evd(w, v)
                    if w is not None:
                        grp_handle = out_handle.create_group(g)
                        grp_handle.create_dataset('tissues', data = np.array(_tissues).astype('S'))
                        grp_handle.create_dataset('w', data = w)
                        grp_handle.create_dataset('v', data = v)
        sys.exit(0)  # this mode stops here
                    
            
        
        

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

