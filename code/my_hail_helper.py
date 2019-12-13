import pandas as pd
import numpy as np
import hail as hl
def names_to_list(names):
    return names.split(',')
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
    pheno = pd.read_csv(csv_path, dtype = type_dic, usecols = list(type_dic.keys())).rename(columns = {indiv_id : 's'})
    covar_table = pheno[covar_list + ['s']]
    trait_table = pheno[pheno_list + ['s']]
    return covar_table, trait_table
def read_indiv_list(file_path):
    # read in the file from file_path
    # and return the first columns as a list of string 
    indiv_list = pd.read_csv(file_path, sep = '\s+', dtype = str)
    return indiv_list.iloc[:,0].to_list()
def subset_by_col(df, colname, target_list):
    return df[df[colname].isin(target_list)]
def df_to_ht(df, key, repartition = 40):
    df = hl.Table.from_pandas(df, key = key)
    df = df.repartition(repartition)
    return df

