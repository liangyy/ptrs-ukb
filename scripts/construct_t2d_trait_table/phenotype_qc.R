args = commandArgs(trailingOnly = TRUE)
library(dplyr)
library(pander)
panderOptions('table.split.table', Inf)
source('../../code/rlib_doc.R')

# load data after population QC
dat = read.csv(args[1])  #  readRDS(args[1])

# code for t2d is 1223: see http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=6
t2d_code = 1223

# first of all, 
# let's see how many instances and how many arraies each trait has.
trait_meta = list()
cols_relavant = colnames(dat)
cols_relavant = cols_relavant[!is.na(stringr::str_match(cols_relavant, '_x_'))]
parsed_col = as.data.frame(do.call(rbind, strsplit(cols_relavant, '_x_')))
colnames(parsed_col) = c('trait', 'instance', 'array')
parsed_col = parsed_col %>% mutate(instance = stringr::str_remove(instance, 'instance_'), array = stringr::str_remove(array, 'array_'))
parsed_col = parsed_col %>% filter(trait != 'ethnicity')  # remove ethnicity since this part is taken care of in population_assignment.Rmd
parsed_col %>% group_by(trait) %>% summarize(ninstance = length(unique(instance))) %>% pander(caption = 'number of instances')
trait_instance_with_multi_array = parsed_col %>% group_by(trait, instance) %>% summarize(number_of_array = length(unique(array))) %>% ungroup() %>% filter(number_of_array > 1) 
message('trait/instance pair with more than one array')
trait_instance_with_multi_array %>% pander

# for T2D, we look into the traits start with self_reported_t2d
# and for each individual, we set to 1 if there is at least one entry = 1223
# otherwise, it is 0
t2d_cols_ind = !is.na(stringr::str_match(colnames(dat), '^self_reported_t2d'))
t2d = apply(dat[, t2d_cols_ind], 1, function(x) {
  (sum(x == t2d_code, na.rm = T) > 0) * 1
})


# for hba1c, we do the following.
# we proceed to aggregate across multiple instances.
# And here we use the first non-NA value as the value of the phenotype for that individual.
trait = 'hba1c'
for(i in trait) {
  this = parsed_col %>% filter(trait == i)
  suffix = '_x_array_0'
  if(i %in% trait_instance_with_multi_array$trait) {
    suffix = '_x_array_agg'
  }
  prefix = paste0(i, '_x_', 'instance_')
  cols = paste0(prefix, unique(this$instance), suffix)
  newcol = data.frame(x = aggregate_instances(dat[, cols], first_non_na))
  hba1c = newcol$x
}

out = data.frame(eid = dat$eid, t2d = t2d, hba1c = hba1c)

message('filter out missing values and HbA1c conc. > 100')
out = out[ rowSums(is.na(out)) == 0, ]
out = out[ out$hba1c <= 100, ]

ptrs_data = read.csv('../../output/query_phenotypes_cleaned_up.csv')
out = inner_join(ptrs_data, out, by = 'eid')

message('quick summary and save results')
out %>% summary
write.csv(out, args[2], quote = F, row.names = F)

message('save the eid by ancestry group')
for(kk in unique(out$meaning)) {
  tmp = out$eid[out$meaning == kk]
  tmp = data.frame(FID = tmp, IID = tmp)
  write.table(tmp, paste0(args[3], kk, '.txt'), row.names = F, col.names = T, quote = F, sep = '\t')
}
