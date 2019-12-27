import pandas as pd
import numpy as np
# import hail as hl
def names_to_list(names):
    if names is None:
        return []
    elif isinstance(names, str):
        return names.split(',')
    else:
        return []
def get_dtype_dic(int_names, str_names, all_names):
    out_dic = { a : np.float for a in all_names }
    for i in int_names:
        out_dic[i] = np.int
    for i in str_names:
        out_dic[i] = str
    return out_dic
def read_and_split_phenotype_csv(csv_path, pheno_names, covar_names, indiv_id, int_names, str_names):
    # read in the phenotype table prepared from ukbREST output (in CSV format)
    # split and return two tables
    # 1. covariate table
    # 2. trait table
    covar_list = names_to_list(covar_names)
    pheno_list = names_to_list(pheno_names)
    all_colnames = covar_list + pheno_list + names_to_list(indiv_id)
    type_dic = get_dtype_dic(names_to_list(int_names), names_to_list(str_names), all_colnames)
    pheno = pd.read_csv(csv_path, dtype = type_dic, usecols = list(type_dic.keys())).rename(columns = {indiv_id : 'eid'})
    covar_table = pheno[covar_list + ['eid']]
    trait_table = pheno[pheno_list + ['eid']]
    return covar_table, trait_table
def read_indiv_list(file_path):
    # read in the file from file_path
    # and return the first columns as a list of string 
    indiv_list = pd.read_csv(file_path, sep = '\s+', dtype = str)
    return indiv_list.iloc[:,0].to_list()
def subset_by_col(df, colname, target_list):
    return df[df[colname].isin(target_list)]
def df_to_ht(df, key, repartition = 40):
    import hail as hl
    df = hl.Table.from_pandas(df, key = key)
    df = df.repartition(repartition)
    return df
def gwas_formater_from_neale_lab(gwas_out, outer_i, inner_j):
    import hail as hl
    # format to GWAS sum stats format from Neale's lab
    # code source: https://github.com/Nealelab/UK_Biobank_GWAS/blob/95ac260a5d4cf9c40effff13fa33fb95ed825e2a/0.2/export_results.biomarkers.py
    i = outer_i
    j = inner_j
    ht_export = gwas_out.annotate(
        n_complete_samples = gwas_out['n'][i],
        AC = gwas_out['sum_x'][i],
        ytx = gwas_out['y_transpose_x'][i][j],
        beta = gwas_out['beta'][i][j],
        se = gwas_out['standard_error'][i][j],
        tstat = gwas_out['t_stat'][i][j],
        pval = gwas_out['p_value'][i][j],
        rsid = gwas_out['rsid'],
        ref = gwas_out['alleles'][0],
        alt = gwas_out['alleles'][1]
    )
    ht_export = ht_export.annotate(
        AF = ht_export['AC'] / (2 * ht_export['n_complete_samples'])
    )
    ht_export = ht_export.annotate(
        minor_AF = hl.cond(ht_export['AF'] <= 0.5, ht_export['AF'], 1.0 - ht_export['AF']),
        minor_allele = hl.cond(ht_export['AF'] <= 0.5, ht_export['alleles'][1], ht_export['alleles'][0])
    )
    ht_export = ht_export.annotate(
        low_confidence_variant = ht_export['minor_AF'] < 0.001
    )
    ht_export = ht_export.select(
        'minor_allele',
        'minor_AF',
        'low_confidence_variant',
        'n_complete_samples',
        'AC',
        'ytx',
        'beta',
        'se',
        'tstat',
        'pval',
        'rsid',
        'ref',
        'alt'
    )
    return ht_export
