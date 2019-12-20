library(optparse)

option_list <- list(
    make_option(c("-i", "--pheno_table"), type="character", default=NULL,
                help="input phenotype table (in ubkREST CSV format)",
                metavar="character"),
    make_option(c("-p", "--indiv_col"), type="character", default=NULL,
                help="column name of individual id",
                metavar="character"),
    make_option(c("-r", "--pop_col"), type="character", default=NULL,
                help="column name of population assignment",
                metavar="character"),
    make_option(c("-q", "--pop"), type="character", default=NULL,
                help="population we want to select from",
                metavar="character"),
    make_option(c("-n", "--nindiv"), type="numeric", default=NULL,
                help="number of individuals desired",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="output with new additional column indicating effect size and standard error converted from OR and p-value",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

set.seed(2019)
df = read.csv(opt$pheno_table)
df = df[df[, opt$pop_col] == opt$pop, ]
idx = sample(1 : nrow(df), size = opt$nindiv, replace = FALSE)
out = df[idx, ]
write.table(out[, opt$indiv_col], opt$output, col = F, row = F, quo = F)
