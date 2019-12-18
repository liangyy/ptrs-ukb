from scipy.stats import norm
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
    return np.apply_along_axis(standardize_vec, 1, mat)
def standardize_vec(vec, offset = 1):
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