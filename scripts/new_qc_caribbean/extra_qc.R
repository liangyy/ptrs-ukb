library(dplyr)
outdir = '../../output'
df = data.table::fread('../../output/new_query_phenotypes_cleaned_up.csv', sep = ',', data.table = F)
d1 = df %>% filter(meaning %in% c('African'))
cores = d1 %>% filter(pc3 <= 0, pc3 >= -12.5, pc4 >= -12.5, pc4 <= 12.5)
cores2 = cores %>% filter(pc1 > 350)
write.table(cores2 %>% select(eid), paste0(outdir, '/CarAfr.txt'), col = T, row = F, quo = F)