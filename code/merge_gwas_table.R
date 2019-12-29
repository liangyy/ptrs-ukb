library(optparse)

option_list <- list(
    make_option(c("-i", "--input_prefix"), type="character", default=NULL,
                help="Prefix of input GWAS table",
                metavar="character"),
    make_option(c("-n", "--input_suffix"), type="character", default=NULL,
                help="Suffix of input GWAS table",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output merged table",
                metavar="character"),
    make_option(c("-t", "--tags"), type="character", default=NULL,
                help="The list of target tags separated by ','",
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
  df = fread(cmd, data.table = F, ...)
  df
}

tags = strsplit(opt$tags, ',')[[1]]
beta_cols = c()
pval_cols = c()
if(length(tags) > 0) {
  filename = paste0(opt$input_prefix, tags[1], opt$input_suffix)
  df = myfread(filename, header = T, sep = '\t')
  df = df[, c('variant', 'beta', 'pval')]
  colnames(df)[2 : 3] = paste0(colnames(df)[2 : 3], '_', tags[1])
  beta_cols = c(beta_cols, colnames(df)[2])
  pval_cols = c(pval_cols, colnames(df)[3])
} else {
  message('No tag. Exit.')
  quit()
}

if(length(tags) > 1) {
  for(t in 2 : length(tags)) {
    filename = paste0(opt$input_prefix, tags[t], opt$input_suffix)
    tmp = myfread(filename, header = T, sep = '\t')
    tmp = tmp[, c('variant', 'beta', 'pval')]
    colnames(tmp)[2 : 3] = paste0(colnames(tmp)[2 : 3], '_', tags[t])
    beta_cols = c(beta_cols, colnames(tmp)[2])
    pval_cols = c(pval_cols, colnames(tmp)[3])
    df = full_join(df, tmp, by = 'variant')
  }
}

for(b in beta_cols) {
  df[is.na(df[, b]), b] = 0
}

for(p in pval_cols) {
  df[is.na(df[, p]), p] = 1
}

write.table(df, opt$output, col = T, row = F, quo = F, sep = '\t')
