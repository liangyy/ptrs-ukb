library(data.table)
library(dplyr)
tmp1 = fread('output/query_first_attempt_with_qc.csv', header = T, sep = ',', data.table = F)
tmp2 = fread('output/query_phenotypes_cleaned_up.csv', header = T, sep = ',', data.table = F)
cols = colnames(tmp1)[c(1, 15:32)]
tmp1 = tmp1[, cols]
tmp2 = tmp2[, cols]
# sum(tmp1 != tmp2)
both_eid = intersect(tmp1$eid, tmp2$eid)
tmp1_not_in_tmp2 = tmp1$eid[which(!tmp1$eid %in% tmp2$eid)]
tmp2_not_in_tmp1 = tmp2$eid[which(!tmp2$eid %in% tmp1$eid)]
tmp1 = tmp1[match(both_eid, tmp1$eid), ]
tmp2 = tmp2[match(both_eid, tmp2$eid), ]
sum(tmp1 != tmp2)

f1 = readRDS('output/query_first_attempt_with_population_qc.rds')
f2 = readRDS('output/query_phenotypes_with_population_qc.rds')
f1$dat_ll %>% filter(eid %in% tmp2_not_in_tmp1)
f2$dat_ll %>% filter(eid %in% tmp1_not_in_tmp2)

a1 = fread('~/Desktop/ukbiobank_query_2/query_2_out.csv', header = T, sep = ',', data.table = F)
a2 = fread('output/query_phenotypes_output.csv', header = T, sep = ',', data.table = F)
a1$eid[which(!a1$eid %in% a2$eid)]
a2$eid[which(!a2$eid %in% a1$eid)]

a1 %>% filter(eid %in% tmp2_not_in_tmp1)
a2 %>% filter(eid %in% tmp1_not_in_tmp2)