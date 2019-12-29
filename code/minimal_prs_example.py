import hail as hl
hl.init()

def read_gwas_table(gwas_table_tsv, type_dic):
    gwas_tsv = hl.import_table(gwas_table_tsv, types = type_dic)
    gwas_tsv = gwas_tsv.annotate(v = hl.parse_variant(gwas_tsv.variant))
    gwas_tsv = gwas_tsv.key_by(gwas_tsv.v.locus, gwas_tsv.v.alleles)
    gwas_tsv = gwas_tsv.drop('v') 
    return gwas_tsv

pval_thresholds = [5e-8, 1e-7, 1e-6, 1e-5]
prefix = 'output_prefix'
subset = 1
mt_subset = hl.read_matrix_table(f'path_to_genotype_{subset}.mt')
for gwas in list(myinputs[subset]['GWASs']):
    gwas_file = myinputs[subset]['GWASs'][gwas]['sum_stat']
    clump_file = myinputs[subset]['GWASs'][gwas]['ld_clump']
    gwas_tsv = read_gwas_table(gwas_file, type_dic = {'beta' : hl.tfloat, 'pval' : hl.tfloat})
    annot_gwas = {
        f'beta_{subset}_x_{gwas}': gwas_tsv[mt_subset.locus, mt_subset.alleles].beta,
        f'pval_{subset}_x_{gwas}': gwas_tsv[mt_subset.locus, mt_subset.alleles].pval
    }
    mt_subset = mt_subset.annotate_rows(**annot_gwas)
    prs = {
        f'pval_thres_{subset}_x_{gwas}_x_{i}' : hl.agg.sum(mt_subset[f'beta_{subset}_x_{gwas}'] * mt_subset.dosage * hl.int(mt_subset[f'pval_{subset}_x_{gwas}'] < i)) for i in pval_thresholds
    }
    mt_subset = mt_subset.annotate_cols(**prs)
    
mt_subset.write(f'{prefix}_{subset}.mt'), overwrite = True)
    