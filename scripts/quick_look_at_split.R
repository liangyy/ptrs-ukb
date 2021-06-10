# quick look at valiation with split
library(dplyr)
library(ggplot2)
options(stringsAsFactors = F)
df = read.csv('~/Desktop/tmp/ptrs-tf/from_nucleus/elastic_net_ptrs_gtex_british_updated_indivs_w_split.performance.csv')
df2 = read.csv('~/Desktop/tmp/ptrs-tf/from_nucleus/prs_british_updated_indivs_w_split.performance.csv')
to_remove = c('basophil', 'mchc')
df = df %>% filter(!trait %in% to_remove)
df2 = df2 %>% filter(!trait %in% to_remove)
df2$sample[df2$sample == 'British_validation'] = 'British_valid'
# head(df)
# head(df2)
myfun = function(df) {
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
  # kk %>% ggplot() + geom_boxplot(aes(x = sample, y = r2, color = sample)) + facet_wrap(~trait, scales = 'free')
  
  kk2 = kk %>% filter(grp == 'repeat0')
  kk %>% ungroup()
}
ptrs = myfun(df)
prs = myfun(df2)

ptrs %>% ggplot() + geom_boxplot(aes(x = sample, y = r2, color = sample)) + facet_wrap(~trait, scales = 'free')
prs %>% ggplot() + geom_boxplot(aes(x = sample, y = r2, color = sample)) + facet_wrap(~trait, scales = 'free')
calc_port = function(dd) {
  dd_s = dd %>% group_by(trait, sample, alpha) %>% summarize(r2_mean = mean(r2), r2_sd = sd(r2)) %>% ungroup() 
  ref = dd_s %>% filter(sample == 'British_valid')
  others = dd_s %>% filter(sample != 'British_valid') %>% left_join(ref %>% select(trait, r2_mean, r2_sd), by = 'trait', suffix = c('', '.ref')) %>% mutate(portability = r2_mean / r2_mean.ref)
}

ptrs_s = calc_port(ptrs)
prs_s = calc_port(prs)
mm = rbind(
  ptrs_s %>% mutate(method = 'ptrs'), 
  prs_s %>% mutate(method = 'prs')
)
mm %>% reshape2::dcast(trait + sample ~ method, value.var = 'portability') %>% ggplot() + geom_point(aes(x = prs, y = ptrs)) + facet_grid(.~sample, scales = 'free') + geom_abline(slope = 1, intercept = 0)
mm %>% ggplot() + geom_boxplot(aes(x = sample, y = portability, color = method))

tmp = read.delim2(pipe('pbpaste'))
tmp$partial_R2 = as.numeric(tmp$partial_R2)
tmp2 = tmp %>% left_join(tmp %>% filter(population == 'EUR ref.'), by = 'trait')
tmp2 = tmp2 %>% mutate(portability = partial_R2.x / partial_R2.y)
tmp2 %>% ggplot() + geom_boxplot(aes(x = population.x, y = portability))
