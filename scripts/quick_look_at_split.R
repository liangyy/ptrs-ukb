# quick look at valiation with split
library(dplyr)
library(ggplot2)
options(stringsAsFactors = F)
df = read.csv('~/Desktop/tmp/ptrs-tf/from_nucleus/elastic_net_ptrs_gtex_british_updated_indivs_w_split.performance.csv')
head(df)
# focus on repeat0
repeats = c('repeat0_1', 'repeat0_2')
# df = df %>% filter(split_label %in% repeats)
df = df %>% mutate(grp = unlist(lapply(strsplit(split_label, '_'), function(x) {x[1]})))
tune_grp = 'repeat0_1'
eval_grp = 'repeat0_2'
get_split_r2 = function(clambda, csplit, cr2, tune_grp, eval_grp) {
  # partial_r2, split_label, lambda
  dd = data.frame(lambda = clambda, split_label = csplit, partial_r2 = cr2)
  tmp = dd %>% reshape2::dcast(lambda ~ split_label, value.var = 'partial_r2')
  x1 = tmp[[tune_grp]]
  x1[is.na(x1)] = -Inf
  x2 = tmp[[eval_grp]]
  x2[is.na(x2)] = -Inf
  res = x2[ which.max(x1)]
  res
}
kk = df %>% group_by(trait, sample, alpha, grp) %>% 
  summarize(
    r2 = get_split_r2(lambda, split_label, partial_r2, tune_grp = paste0(grp[1], '_1'), eval_grp = paste0(grp[1], '_2'))
  )
kk %>% ggplot() + geom_boxplot(aes(x = sample, y = r2, color = sample)) + facet_wrap(~trait, scales = 'free')

kk2 = kk %>% filter(grp == 'repeat0')
