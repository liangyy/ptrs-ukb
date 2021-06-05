d1 = data.table::fread('../../output/new_query_phenotypes_cleaned_up.csv', sep = ',', data.table = F)
d2 = data.table::fread('../../output/query_phenotypes_cleaned_up.csv', sep = ',', data.table = F)
dd = inner_join(d1, d2, by = 'eid')
summary(dd$rbc.x - dd$rbc.y)
kk = colnames(d1)
for(k in kk) {
  if(k %in% c('eid', 'meaning')) {
    next
  }
  print(summary(dd[[paste0(k, '.x')]] - dd[[paste0(k, '.y')]]))
}

message('They are the same!')