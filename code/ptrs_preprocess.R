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

library(dplyr)

df = read.csv(opt$spredixcan)

if(opt$mode == 'naive') {
  df = df %>% filter(abs(effect_size) < 5)
  write.table(df$gene, opt$output, quo = F, col = F, row = F)
}
