library(optparse)

option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Input file to extract from",
                metavar="character"),
    make_option(c("-e", "--extract_list"), type="character", default=NULL,
                help="A list of target to extract",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output",
                metavar="character"),
    make_option(c("-x", "--extract_col"), type="character", default=NULL,
                help="Column name in input to look for",
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

dffile = myfread(opt$input, header = T, sep = '\t')
extract_list = myfread(opt$extract_list, header = F, sep = '\t')
sub = dffile[dffile[, opt$extract_col] %in% extract_list$V1, ]
write.table(sub, opt$output, col = T, row = F, sep = '\t', quote = F)
