# generate a list of genes which exist in both MESA CAU and MESA AFHI.

# conda activate predixcan_prediction

mesa_cau_h5 = '/vol/bmd/yanyul/UKB/predicted_expression/predicted_expression.ukb_imp_x_CAU.h5'
mesa_afhi_h5 = '/vol/bmd/yanyul/UKB/predicted_expression/predicted_expression.ukb_imp_x_AFHI.h5'

outlist = '../misc/common_genes_in_mesa_cau_and_afhi.txt'

import h5py

with h5py.File(mesa_cau_h5, 'r') as f:
    gene_cau = f['genes'][:].astype(str)

with h5py.File(mesa_afhi_h5, 'r') as f:
    gene_afhi = f['genes'][:].astype(str)

gene_common = set(gene_cau)
gene_common = gene_common.intersection(set(gene_afhi))
gene_common = list(gene_common)

with open(outlist, 'w') as f:
    for g in gene_common:
        f.write(g +'\n')

    
