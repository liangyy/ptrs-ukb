library(SilverStandardPerformance)
library(dplyr)
data("gene_annotation_gencode_v26_hg38")
data("ld_block_pickrell_eur_b38")
o = SilverStandardPerformance:::annotate_gwas_loci_with_gene(ld_block_pickrell_eur_b38$ld_block %>% mutate(trait = 'placeholder'), gene_annotation_gencode_v26_hg38$gene_annotation)

o = o %>% mutate(ld_block = paste0(chromosome, ':', start, ':', end)) %>% select(gene, ld_block)
write.table(o, '../../output/ld_block_annotation_of_genes.tsv', row = F, col = T, quo = F, sep = '\t')
