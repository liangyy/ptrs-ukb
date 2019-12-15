args = commandArgs(trailingOnly = TRUE)
library(dplyr)
library(pander)
panderOptions('table.split.table', Inf)
source('../../code/rlib_doc.R')

# load data after population QC
dat = readRDS(args[1])
dat = dat$dat_pass_pop_qc

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

# for trait/instance pairs with multiple arrays, 
# It seems that taking the values from multiple arraies are quite consistent so that we aggregate array by taking the average.
for(i in 1 : nrow(trait_instance_with_multi_array)) {
  this = trait_instance_with_multi_array[i, ]
  prefix = paste0(this$trait, '_x_', 'instance_', this$instance, '_x_', 'array_')
  cols = paste0(prefix, 0 : (this$number_of_array - 1))
  newcol = data.frame(x = aggregate_instances(dat[, cols], aggregate_instances))
  colnames(newcol) = paste0(prefix, 'agg')
  dat = cbind(dat, newcol)
  dat[, cols] = NULL
}

# then, we proceed to aggregate across multiple instances.
# And here we use the first non-NA value as the value of the phenotype for that individual.
for(i in unique(parsed_col$trait)) {
  this = parsed_col %>% filter(trait == i)
  suffix = '_x_array_0'
  if(i %in% trait_instance_with_multi_array$trait) {
    suffix = '_x_array_agg'
  }
  prefix = paste0(i, '_x_', 'instance_')
  cols = paste0(prefix, unique(this$instance), suffix)
  newcol = data.frame(x = aggregate_instances(dat[, cols], first_non_na))
  colnames(newcol) = i
  dat = cbind(dat, newcol)
  dat[, cols] = NULL
}
martin_traits = as.character(unique(parsed_col$trait))
summary(dat[, martin_traits])

# to make the analysis easier to work with, let's try to limit to individual with all phenotypes being non-missing.
all_none_missing = rowSums(is.na(dat[, martin_traits])) == 0
dat %>% filter(all_none_missing) %>% group_by(meaning) %>% summarize(nindiv = n()) %>% pander(caption = 'Number of individuals after limiting to ones with all phenotypes non-missing')
dat_cleaned = dat %>% filter(all_none_missing) %>% select(-ethnicity_x_instance_0_x_array_0, -ethnicity_x_instance_1_x_array_0, -ethnicity_x_instance_2_x_array_0)

# show the cleaned-up data (after filtering)
message('quick summary before save')
dat_cleaned %>% summary

# add more covariates:
# 1. age_squared
# 2. age_times_sex
# 3. age_squared_times_sex
message('adding more covariates: age x sex ...')
dat_cleaned = dat_cleaned %>% mutate(age_squared = age_recruitment ^ 2, age_times_sex = age_recruitment * sex, age_squared_times_sex = (age_recruitment ^ 2) * sex)

message('quick summary and save results')
dat_cleaned %>% summary
write.csv(dat_cleaned, args[2], quote = F, row.names = F)
