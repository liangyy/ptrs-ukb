library(optparse)

option_list <- list(
    make_option(c("-p", "--prs"), type="character", default=NULL,
                help="PRS file (TSV.BGZ)",
                metavar="character"),
    make_option(c("-i", "--indiv_lists"), type="character", default=NULL,
                help="list of individual lists separated by ','",
                metavar="character"),
    make_option(c("-n", "--indiv_col"), type="character", default=NULL,
                help="column name of individual ID",
                metavar="character"),
    make_option(c("-e", "--pheno_table"), type="character", default=NULL,
                help="Phenotype/covariate table in CSV format",
                metavar="character"),
    make_option(c("-y", "--pheno_yaml"), type="character", default=NULL,
                help="YAML file specifying columns of phenotypes and covariates, etc",
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
  if(ext == 'gz' | ext == 'bgz') {
    cmd = paste0('zcat < ', filename)
  } else {
    cmd = filename
  }
  df = fread(cmd, data.table = F, ...)
  df
}

parse_trait_and_pthreshod_from_prs = function(str) {
  tmp = lapply(strsplit(str, '_x_'), function(x) {
    data.frame(trait = x[2], pthreshold = x[3])
  })
  tmp = do.call(rbind, tmp)
  tmp %>% mutate(colname = str)
}


# read in population assignment
pop_list = list()
for(pop in strsplit(opt$indiv_lists, ',')[[1]]) {
  pop_name = stringr::str_remove(basename(pop), '.txt')
  tmp = myfread(pop, header = F)
  tmp = as.character(tmp$V1)
  pop_list[[length(pop_list)]] = data.frame(indiv = tmp, population = pop_name)
}
pop_assign = do.call(rbind, pop_list)


# read in prs
prs = myfread(opt$prs, header = T, sep = '\t')
prs_cols = colnames(prs)
prs_cols[which(prs_cols == opt$indiv_col)] = NULL
trait_and_pvals = parse_trait_and_pthreshod_from_prs(prs_cols)
join_by = 'indiv'; names(join_by) = opt$indiv_col
message('Before join nrow(prs) = ', nrow(prs))
prs = inner_join(prs, pop_assign, by = join_by)
message('After join nrow(prs) = ', nrow(prs))


# read phenotype/covariates
pheno = myfread(opt$pheno_table, header = T, sep = ',')
meta_info = yaml::read_yaml(opt$pheno_yaml)
covars = strsplit(meta_info$covar_names, ',')[[1]]
df_covar = pheno[, c(covars, meta_info$indiv_id)]
df_pheno = pheno[, c(unique(trait_and_pvals$trait), meta_info$indiv_id)]
join_by = meta_info$indiv_id
names(join_by) = opt$indiv_col
prs = inner_join(prs, df_covar, by = join_col)
prs = inner_join(prs, df_pheno, by = join_col)


# loop over traits
out = list()
for(i in 1 : nrow(trait_and_pvals)) {
  colname_i = trait_and_pvals$colname[i]
  trait_i = trait_and_pvals$trait[i]
  pval_i = trait_and_pvals$pthreshold[i]
  tmp = prs %>% group_by(population) %>% do(compute_r2(., colname_i, trait_i, covars)
  out[[length(out) + 1]] = tmp %>% mutate(trait = trait_i, pval_cutoff = pval_i)
}
out = do.call(rbind, out)
write.table(out, opt$output, col = T, row = F, quo = F, sep = '\t')

