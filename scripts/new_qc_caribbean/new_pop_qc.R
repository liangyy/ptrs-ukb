args = commandArgs(trailingOnly = TRUE)
library(dplyr)
library(data.table)
library(pander)
panderOptions('table.split.table', Inf)
source('../../code/rlib_doc.R')


# load phenotype table obtained from query
dat = fread(args[1], data.table = F)

# load data coding downloaded from ukb website
# see here https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001
data_coding_for_ethnicity = read.delim2('~/Downloads/coding1001.tsv')

# clean up the coding table
# subset to what we want: african, british, chinese, indian
target_ethnicity_groups = c('Chinese', 'British', 'Indian', 'African', 'Pakistani', 'Bangladeshi', 'Caribbean')
new_meaning_map = data.frame(
  old = c('Chinese', 'British', 'Indian', 'African', 'Pakistani', 'Bangladeshi', 'Caribbean'),
  new = c('Chinese', 'British', 'Indian', 'African', 'Indian', 'Indian', 'Caribbean')
)
data_coding_for_ethnicity = data_coding_for_ethnicity %>% filter(meaning %in% target_ethnicity_groups)

# aggregate ukb ethnicity across instance by using the first non-missing ethinicity 
ethnicity_label = dat %>% select(ethnicity_x_instance_0_x_array_0, ethnicity_x_instance_1_x_array_0, ethnicity_x_instance_2_x_array_0) 
dat$ethnicity_agg = aggregate_instances(ethnicity_label)


# filter out individuals with un-wanted ethnicity
dat = dat %>% filter(!is.na(ethnicity_agg)) %>% filter(ethnicity_agg %in% data_coding_for_ethnicity$coding)
dat = dat %>% inner_join(data_coding_for_ethnicity %>% select(coding, meaning), by = c('ethnicity_agg' = 'coding'))
message('quick look at data')

dat %>% group_by(meaning) %>% summarize(nindiv = n()) %>% pander
# update meaning
tmp = left_join(dat %>% select(meaning), new_meaning_map, by = c('meaning' = 'old')) %>% pull(new)
dat$meaning = tmp
dat %>% group_by(meaning) %>% summarize(nindiv = n()) %>% pander


# extract the first 10 PC for population QC
# and calculate the probability (in log scale) of a individual in each ethnicity cluster (treated as multivariate normal)
vis = dat %>% select(pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10, meaning)
ll_list = list()
for(i in unique(vis$meaning)) {
  sub = vis %>% filter(meaning == i) %>% select(-meaning)
  mean_hat = colMeans(sub)
  cov_hat = cov(sub)
  tmp = data.frame(ll = mvtnorm::dmvnorm(vis %>% select(-meaning), mean = mean_hat, sigma = cov_hat, log = T))
  colnames(tmp) = paste(i, 'dist')
  ll_list[[length(ll_list) + 1]] = tmp
}
df_ll = do.call(cbind, ll_list)
df_ll = df_ll %>% mutate(meaning = vis$meaning, eid = dat$eid)
message('log-probability of individuals in different populations')
df_ll %>% head %>% pander

# perform filtering 
# for individuals in each ethnicity group:
# if the log probability < -50 in the corresponding ethnicity group,
# the individual is removed
ll_cutoff = -50
eid_pass_cutoff = c()
for(pop in unique(df_ll$meaning)) {
  dist = paste(pop, 'dist')
  sub = df_ll %>% filter(meaning == pop)
  eid_pass_cutoff = c(eid_pass_cutoff, sub$eid[sub[, dist] > ll_cutoff])
}
dat_pass_cutoff = dat %>% filter(eid %in% eid_pass_cutoff)

# save the data passing population QC
saveRDS(list(dat_pass_pop_qc = dat_pass_cutoff, dat_ll = df_ll), args[2])