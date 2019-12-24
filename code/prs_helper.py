import hail as hl

def read_gwas_table_with_varlist(gwas_table_tsv, varlist, type_dic):
    # varlist can be one file or a list of files
    # type dic is for gwas_table_tsv
    gwas_tsv = hl.import_table(gwas_table_tsv, key = ['rsid'], types = type_dic)
    clump_snp = hl.import_table(varlist, key = ['f0'], no_header = True)
    gwas_tsv = gwas_tsv.filter(hl.is_defined(clump_snp[gwas_tsv.rsid]))
    k = hl.parse_variant(gwas_tsv.variant)
    gwas_tsv = gwas_tsv.annotate(**k)
    gwas_tsv = gwas_tsv.key_by(gwas_tsv.locus, gwas_tsv.alleles)
    gwas_tsv = gwas_tsv.repartition(40)
    gwas_tsv = gwas_tsv.cache()
    return gwas_tsv