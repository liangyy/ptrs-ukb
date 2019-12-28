import hail as hl
import os

def read_gwas_table_with_varlist(gwas_table_tsv, varlist, type_dic, checkpoint_path, gwas_ht = None, no_return = False):
    # varlist can be one file or a list of files
    # type dic is for gwas_table_tsv
    if gwas_ht is None:
        gwas_tsv = hl.import_table(gwas_table_tsv, key = ['rsid'], types = type_dic)
    else:
        phenotypes = gwas_ht['phenotypes'].collect()[0]
        i, j = get_index_in_nested_list(phenotypes, gwas_table_tsv)
        gwas_tsv = gwas_ht.annotate(
            beta = gwas_ht['beta'][i][j],
            pval = gwas_ht['p_value'][i][j]
        )
        gwas_tsv = gwas_tsv.select(
            'beta',
            'pval',
            'locus',
            'alleles'
        )
        
    clump_snp = hl.import_table(varlist, key = ['f0'], no_header = True)
    gwas_tsv = gwas_tsv.filter(hl.is_defined(clump_snp[gwas_tsv.rsid]))
    if gwas_ht is None:
        k = hl.parse_variant(gwas_tsv.variant)
        gwas_tsv = gwas_tsv.annotate(**k)
    gwas_tsv = gwas_tsv.key_by(gwas_tsv.locus, gwas_tsv.alleles)
    # gwas_tsv = gwas_tsv.repartition(40)
    # gwas_tsv = gwas_tsv.cache()
    if no_return is False:
        gwas_tsv = gwas_tsv.checkpoint(checkpoint_path, overwrite = True)
        return gwas_tsv
    else:
        gwas_tsv.write(checkpoint_path, overwrite = True)
    

def get_index_in_nested_list(inlist, target):
    outer_i = None
    inner_j = None
    for i in range(len(inlist)):
        try:
            outer_i = i
            inner_j = inlist[i].index(target)
            break
        except:
            pass
    # print(outer_i, inner_j)
    return outer_i, inner_j

def remove_ht(filepath):
    os.system("rm -rf {}".format(filepath))
    
    
def read_gwas_table(gwas_table_tsv, type_dic):
    gwas_tsv = hl.import_table(gwas_table_tsv, types = type_dic)
    gwas_tsv = gwas_tsv.annotate(v = hl.parse_variant(gwas_tsv.variant))
    gwas_tsv = gwas_tsv.key_by(gwas_tsv.v.locus, gwas_tsv.v.alleles)
    gwas_tsv = gwas_tsv.drop('v') 
    return gwas_tsv
    