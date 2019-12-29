library(optparse)

option_list <- list(
    make_option(c("-p", "--input_prefix"), type="character", default=NULL,
                help="prefix of input",
                metavar="character"),
    make_option(c("-s", "--input_suffix"), type="character", default=NULL,
                help="suffix of input",
                metavar="character"),
    make_option(c("-n", "--indiv_col"), type="character", default=NULL,
                help="column name of individual ID",
                metavar="character"),
    make_option(c("-l", "--populations"), type="character", default=NULL,
                help="populations separated by :",
                metavar="character"),
    make_option(c("-e", "--pheno_table"), type="character", default=NULL,
                help="Phenotype/covariate table in CSV format",
                metavar="character"),
    make_option(c("-y", "--pheno_yaml"), type="character", default=NULL,
                help="YAML file specifying columns of phenotypes and covariates, etc",
                metavar="character"),
    make_option(c("-t", "--trait_col"), type="character", default=NULL,
                help="column of trait",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output R2",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(dplyr)
library(data.table)
options(stringsAsFactors = FALSE)
source('../../code/rlib_doc.R')


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

pheno = myfread(opt$pheno_table, header = T, sep = ',')
meta_info = yaml::read_yaml(opt$pheno_yaml)
covars = strsplit(meta_info$covar_names, ',')[[1]]
df_covar = pheno[, c(covars, meta_info$indiv_id)]
df_pheno = pheno[, c(opt$trait_col, meta_info$indiv_id)]


pops = strsplit(opt$populations, ':')[[1]]
ptrs = list()
for(p in pops) {
  filename = paste0(opt$input_prefix, p, opt$input_suffix)
  tmp = myfread(filename, header = T, sep = '\t')
  ptrs[[length(ptrs) + 1]] = tmp %>% mutate(population = p)
}
ptrs = do.call(rbind, ptrs)
ptrs_cols = colnames(ptrs)
ptrs_cols = ptrs_cols[!ptrs_cols %in% c(opt$indiv_col, 'population')]

join_col = meta_info$indiv_id
names(join_col) = opt$indiv_col
ptrs = inner_join(ptrs, df_covar, by = join_col)
ptrs = inner_join(ptrs, df_pheno, by = join_col)

out = list()
for(ptrs_k in ptrs_cols) {
  tmp = ptrs %>% group_by(population) %>% do(report_r2(., opt$trait_col, ptrs_k, covars))
  out[[length(out) + 1]] = tmp %>% mutate(ptrs_col = ptrs_k)
}
out = do.call(rbind, out)
write.table(out, opt$output, col = T, row = F, quo = F, sep = '\t')

