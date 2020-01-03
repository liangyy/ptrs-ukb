library(optparse)

option_list <- list(
    make_option(c("-s", "--spredixcan"), type="character", default=NULL,
                help="S-PrediXcan output",
                metavar="character"),
    make_option(c("-p", "--pred_expr"), type="character", default=NULL,
                help="gene by individual matrix (TSV.GZ)",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output PTRS",
                metavar="character"),
    make_option(c("-e", "--gene_list"), type="character", default=NULL,
                help="Gene list",
                metavar="character"),
    make_option(c("-n", "--gene_list_pval"), type="character", default='No',
                help="If want to use pval contained in gene list, set to 'Yes'",
                metavar="character"),
    make_option(c("-v", "--pval_cutoffs"), type="character", default='1e-7,0.05',
                help="pvalue cutoffs separated by ','",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(dplyr)
library(data.table)
options(stringsAsFactors = FALSE)


myfread = function(filename, ...) {
  ext = tools::file_ext(filename)
  if(ext == 'gz') {
    cmd = paste0('zcat < ', filename)
  } else {
    cmd = filename
  }
  df = fread(cmd, ...)
  df
}

spredixcan = read.csv(opt$spredixcan)
pred_expr = myfread(opt$pred_expr, header = T, sep = '\t', data.table = F)
gene_list = read.table(opt$gene_list, header = F)$V1
pvals = as.numeric(unlist(strsplit(opt$pval_cutoffs, ',')[[1]]))

gene_pool = intersect(pred_expr$gene, spredixcan$gene)
gene_pool = intersect(gene_list, gene_pool)

pred_expr = pred_expr[pred_expr$gene %in% gene_pool, ]
pred_expr_mat = pred_expr[, which(colnames(pred_expr) != 'gene')]
indiv_list = colnames(pred_expr_mat)
gene_list = pred_expr$gene
pred_expr_mat = t(as.matrix(pred_expr_mat))

spredixcan_cleaned = spredixcan %>% filter(gene %in% gene_pool)
spredixcan_cleaned = spredixcan_cleaned[match(gene_list, spredixcan_cleaned$gene), ]

if(opt$gene_list_pval == 'Yes') {
  gene_list_df = read.table(opt$gene_list, header = F, sep ='\t')
  colnames(gene_list_df) = c('gene', 'pvalue')
  spredixcan_cleaned = inner_join(spredixcan_cleaned %>% select(-pvalue), gene_list_df, by = 'gene')
}

gene_beta = spredixcan_cleaned$effect_size
gene_pval = spredixcan_cleaned$pvalue

# message(dim(pred_expr_mat))
# message(length(gene_beta))

out = list()
out[[length(out) + 1]] = data.frame(x = indiv_list)
for(i in pvals) {
  ptrs = pred_expr_mat %*% (gene_beta * (gene_pval < i))
  out[[length(out) + 1]] = data.frame(x = ptrs)
}
out = do.call(cbind, out)
colnames(out) = c('indiv', paste0('pval_', pvals))
gz1 <- gzfile(opt$output, "w")
write.table(out, gz1, col = T, row = F, quo = F, sep = '\t')
close(gz1)
