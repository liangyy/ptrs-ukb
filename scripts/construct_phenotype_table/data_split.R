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

# Let's keep 5000 British as validation set and 5000 British as test set. 
# For each split file, we keep a three lists called: `c.N.txt`, `british-validation.N.txt`, and `british-test.N.txt` along with 3 lists for other 3 populations. 
# do the splitting here
set.seed(2019)
ntest = 5000
nvalid = 5000
outdir = args[2]
if(!dir.exists(outdir)) {
  dir.create(outdir)
}
for(i in 1 : nsplits) {
  df_b = df %>% filter(meaning == 'British')
  df_b = data.frame(FID = df_b$eid, IID = df_b$eid)
  testvalididx = sample(1 : nrow(df_b), size = ntest + nvalid, replace = F)
  train = df_b[! (1 : nrow(df_b)) %in% testvalididx, ]
  sub = df_b[testvalididx, ]
  testidx = sample(1 : (ntest + nvalid), size = ntest, replace = F)
  test = sub[ (1 : nrow(sub)) %in% testidx, ]
  valid = sub[! (1 : nrow(sub)) %in% testidx, ]
  write.table(train, paste0(outdir, '/', 'British-training-', i, '.txt'), col = T, row = F, quo = F)
  write.table(test, paste0(outdir, '/', 'British-test-', i, '.txt'), col = T, row = F, quo = F)
  write.table(valid, paste0(outdir, '/', 'British-validation-', i, '.txt'), col = T, row = F, quo = F)
}
for(p in unique(df$meaning)) {
  if(p == 'British') {
    next
  }
  sub = df %>% filter(meaning == p)
  sub = data.frame(FID = sub$eid, IID = sub$eid)
  write.table(sub, paste0(outdir, '/', p, '.txt'), col = T, row = F, quo = F)
}

