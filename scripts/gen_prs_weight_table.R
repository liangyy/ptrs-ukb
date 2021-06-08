# Generate PRS weight table
# format:
# rsid effect_allele weight1 weights2 ...
library(optparse)

option_list <- list(
  make_option(c("-d", "--data_yaml"), type="character", default=NULL,
              help="data yaml: 
              [trait-name]: 
                [gwas-table]: 'path'
                [snp-pool]: 'path'
              ",
              metavar="character"),
  # make_option(c("-p", "--pval_cutoffs"), type="character", default=NULL,
  #             help="List of p-value cutoffs, separated by ,",
  #             metavar="character"),
  make_option(c("-l", "--log_level"), type="integer", default=9,
              help="Logging level",
              metavar="character"),
  make_option(c("-n", "--snpid"), type="character", default=NULL,
              help="column name of SNP ID",
              metavar="character"),
  make_option(c("-v", "--pval"), type="character", default=NULL,
              help="column name of p-values",
              metavar="character"),
  make_option(c("-e", "--effect_size"), type="character", default=NULL,
              help="column name of effect size",
              metavar="character"),
  make_option(c("-f", "--effect_allele"), type="character", default=NULL,
              help="column name of effect allele",
              metavar="character"),
  make_option(c("-o", "--output_prefix"), type="character", default=NULL,
              help="Prefix of output",
              metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser) 

library(data.table)
library(dplyr)

# logging config
logging::basicConfig(level = opt$log_level)

data_list = yaml::read_yaml(opt$data_yaml)
# pval_cutoffs = as.numeric(strsplit(opt$pval_cutoffs, ',')[[1]])

for(trait in names(data_list)) {
  
  gwas_tsv = data_list[[trait]]$sum_stat
  snp_pool = data_list[[trait]]$ld_clump
  
  logging::loginfo(glue::glue('{trait}: Loading GWAS table'))
  df_gwas = fread(gwas_tsv, sep = '\t', data.table = F)
  
  logging::loginfo(glue::glue('{trait}: Loading SNP Pool'))
  snp_pool = fread(snp_pool, header = F, data.table = F)$V1
  
  logging::loginfo(
    glue::glue('{trait}: There are {nrow(df_gwas)} SNPs in GWAS.')
  )
  
  logging::loginfo(
    glue::glue('{trait}: There are {length(snp_pool)} SNPs in SNP pool')
  )
  
  df_in = df_gwas[ df_gwas[[opt$snpid]] %in% snp_pool, ]
  logging::loginfo(
    glue::glue('{trait}: There are {nrow(df_in)} SNPs in both')
  )
  df_in = df_in[, c(opt$snpid, opt$pval, opt$effect_allele, opt$effect_size)] 
  colnames(df_in) = c('snpid', 'pval', 'effect_allele', 'effect_size')
  # df_weights = NULL
  # for(pp in sort(pval_cutoffs, decreasing = T)) {
  #   df_in = df_in[ df_in[[opt$pval]] <= pp, ]
  #   logging::loginfo(
  #     glue::glue('{trait}: There are {nrow(df_in)} SNPs at p-value <= {pp}')
  #   )
  #   tmp = df_in[, c(opt$snpid, opt$effect_allele, opt$effect_size)]
  #   colnames(tmp) = c('snpid', 'effect_allele', glue::glue('pval_thres_x_{trait}_x_{pp}'))
  #   if(is.null(df_weights)) {
  #     df_weights = tmp
  #   } else {
  #     df_weights = left_join(df_weights, tmp, by = c('snpid', 'effect_allele'))
  #   }
  # }
  # df_weights[is.na(df_weights)] = 0
  
  logging::loginfo(
    glue::glue('{trait}: Saving')
  )
  fn = glue::glue('{opt$output_prefix}.{trait}.prs_weights.tsv')
  write.table(df_in, fn, row = F, col = T, quo = F, sep = '\t')
}



