from scipy.stats import norm
import numpy as np
import pandas as pd
def extract_cols(h5mat, col_names, target_names):
    col_names = np.array(col_names)
    target_names = np.array(target_names)
    target_idx = np.where(np.isin(col_names, target_names))[0].flatten().tolist()
    subcols = col_names[target_idx]
    submat = h5mat[:, target_idx]
    return submat, subcols
def remove_constant_row(mat, row_names):
    row_names = np.array(row_names)
    std = np.std(mat, axis = (1))
    return mat[std != 0, :], row_names[std != 0]
def inv_norm_row(mat):
    return np.apply_along_axis(inv_norm_vec, 1, mat)
def inv_norm_vec(vec, offset = 1):
    rank = myrank(vec)
    return norm.ppf(rank / (len(rank) + offset), loc = 0, scale = 1)
def myrank(vec):
    argsort = np.argsort(vec)
    ranks = np.empty_like(argsort)
    ranks[argsort] = np.arange(len(vec))
    return ranks + 1  # rank starts from 1
def format_to_gcta_grm(mat):
    nobs = mat.shape[0]
    nindiv = mat.shape[1]
    grm = np.dot(mat.transpose(), mat) / nobs
#     row_index = [ i for j in range(nindiv - i) for i in range(nindiv) ]  # 0-base
#     col_index = [ j for j in range(i, nindiv) for i in range(nindiv) ]  # 0-base
    indices = np.tril_indices(nindiv)
    grm_values = grm[indices]
    nobs_vec = np.ones(grm_values.shape[0]) * nobs
    df = pd.DataFrame({'idx1': indices[0] + 1, 'idx2': indices[1] + 1, 'nobs': nobs_vec.astype(int), 'gr': grm_values})
    return df
def get_header(filename):
    f = open(filename, 'r')
    samples = f.readline().strip().split('\t')
    f.close()
    return samples
def read_pred_expr_as_mt(filename):
    # the file should be TSV
    # it should have a column named 'gene'
    # all other columns are considered as samples
    # it returns hail MatrixTable with columns as samples (key 's') and rows as gene (key 'gene')
    # and entry 'pred_expr'
    samples = get_header(filename)
    samples.pop(samples.index('gene'))  # we remove 'gene' here since we only want samples
    type_dic = {
        k : hl.tfloat for k in samples
    }
    type_dic['gene'] = hl.tstr
    pred_expr = hl.import_table(filename, types = type_dic)
    pred_expr = (pred_expr
        .key_by('gene')
        .to_matrix_table_row_major(columns = samples, entry_field_name = 'pred_expr', col_field_name = 's')
    )
    return pred_expr
def read_tsv(filename, sample_col):
    # the file should be TSV
    # it should have one column named sample_col for sample ID
    # all other columns are considered as float!
    # it returns hail Table with row as sample (key 's') 
    float_cols = get_header(filename)
    float_cols.pop(float_cols.index(sample_col))  # we remove the column for sample ID here
    type_dic = {
        k : hl.tfloat for k in float_cols
    }
    type_dic[sample_col] = hl.tstr
    ht = hl.import_table(filename, types = type_dic)
    ht = ht.annotate(s = ht[sample_col])
    ht = ht.key_by('s')
    ht = ht.drop(sample_col)
    ht = ht.repartition(ht.n_partitions())
    return ht
def struct_to_np_array(struct, exclude_cols):
    m = len(struct[0].collect())
    colnames = list(struct.keys())
    n = np.sum(np.logical_not(np.isin(colnames, exclude_cols)))
    out_array = np.empty([m, n])
    counter = 0
    for i in range(len(colnames)):
        if colnames[i] not in exclude_cols:
            out_array[:, i] = struct[i].collect()
    return out_array
def tsv_to_pd_df(filename, indiv_col):
    return  pd.read_csv(filename, header = 0, sep = '\t', dtype = {indiv_col: str})
def standardize_row(mat):
    return np.apply_along_axis(standardize_vec, 1, mat)
def standardize_vec(vec):
    return _divide(vec - np.mean(vec), np.std(vec))
def _divide(a, b):
    return np.divide(a, b, out = np.zeros_like(a), where = (b! = 0))
def truncate_evd(w, v, lambda_max_over_lambda = 1 / 30):
    w_max = np.max(w)
    w_scaled = w / w_max
    idx = np.where(w_scaled > lambda_max_over_lambda)[0]
    if idx.size == 0:
        w_keep = None
        v_keep = None
    else:
        idx = np.min(idx)
        w_keep = w[idx:]
        v_keep = v[:, idx:]
    return w_keep, v_keep