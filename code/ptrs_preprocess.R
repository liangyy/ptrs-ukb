library(optparse)

option_list <- list(
    make_option(c("-s", "--spredixcan"), type="character", default=NULL,
                help="S-PrediXcan output",
                metavar="character"),
    make_option(c("-m", "--mode"), type="character", default='naive',
                help="Preprocessing mode (only naive is implemented)",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output gene list",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

source('../../code/rlib_doc.R')

library(dplyr)
options(stringsAsFactors = F)

df = read.csv(opt$spredixcan)

if(opt$mode == 'naive') {
  df = df %>% filter(abs(effect_size) < 5)
  write.table(df$gene, opt$output, quo = F, col = F, row = F)
} else if(opt$mode == 'ldblock') {
  annot = read.table('../../output/ld_block_annotation_of_genes.tsv', header = T, sep = '\t')
  df = df %>% mutate(gene_no_dot = trim_dot(gene))
  df = inner_join(df, annot, by = c('gene_no_dot' = 'gene'))
  df_top = df %>% group_by(ld_block) %>% summarize(top_significance = max(abs(zscore))) %>% ungroup()
  df = df %>% inner_join(df_top, by = 'ld_block') %>% filter(abs(zscore) == top_significance)
  df = df %>% filter(abs(effect_size) < 5)
  write.table(df$gene, opt$output, quo = F, col = F, row = F)
}
