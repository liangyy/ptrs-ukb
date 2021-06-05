args = commandArgs(trailingOnly = TRUE)

library(dplyr)
library(pander)
library(data.table)
panderOptions('table.split.table', Inf)
source('../../code/rlib_doc.R')

# load the phenotype table
df = fread(args[1], sep = ',', header = T, data.table = F)
df_paper = read.csv('../../external_data/martin_et_al_2019ng_table_s6.csv')
traits = tolower(df_paper$Trait)
nsplits = length(traits)

# just to see the number of individuals in each population again
message('nindiv by population')
df %>% group_by(meaning) %>% summarize(nindiv = n()) %>% pander


set.seed(2019)
outdir = args[2]
if(!dir.exists(outdir)) {
  dir.create(outdir)
}
# skip EUR

for(p in unique(df$meaning)) {
  if(p == 'British') {
    next
  }
  sub = df %>% filter(meaning == p)
  sub = data.frame(FID = sub$eid, IID = sub$eid)
  write.table(sub, paste0(outdir, '/', p, '.txt'), col = T, row = F, quo = F)
}

