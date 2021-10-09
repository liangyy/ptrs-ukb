# PTRS performance 
# Baboon genes vs all genes
library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 20))
shift = 0.1
source('code/rlib_doc.R')
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
df1 = load_perf('~/Desktop/tmp/ptrs-tf/from_nucleus/elastic_net_ptrs_gtex_british_split_british_w_split_lambda.performance.csv') %>% filter(!trait %in% bad_traits)
df1pt = load_perf('~/Desktop/tmp/ptrs-tf/from_nucleus/pt_ptrs_gtex_british_split_british_w_split_lambda.performance.csv') %>% filter(!trait %in% bad_traits)
df2 = load_perf('~/Desktop/tmp/ptrs-tf/from_nucleus/elastic_net_ptrs_gtex_british_split_british_baboon_genes_w_split_lambda.performance.csv') %>% filter(!trait %in% bad_traits)
df2pt = load_perf('~/Desktop/tmp/ptrs-tf/from_nucleus/pt_ptrs_gtex_british_split_british_baboon_genes_w_split_lambda.performance.csv') %>% filter(!trait %in% bad_traits)
df1s = get_test_perf_from_splits(df1)
df1spt = get_test_perf_from_splits(df1pt)
df2s = get_test_perf_from_splits(df2)
df2spt = get_test_perf_from_splits(df2pt)
df1ss = summarize_across_grps(df1s)
df1sspt = summarize_across_grps(df1spt)
df2ss = summarize_across_grps(df2s)
df2sspt = summarize_across_grps(df2spt)

df1ssp = calc_port(df1ss)
df1ssppt = calc_port(df1sspt)
df2ssp = calc_port(df2ss)
df2ssppt = calc_port(df2sspt)


tmp = rbind(
  df1ssp %>% mutate(method = 'PTRS'),
  df1ssppt %>% mutate(method = '(PT) PTRS'),
  df2ssp %>% mutate(method = 'PTRS-b'),
  df2ssppt %>% mutate(method = '(PT) PTRS-b')
) %>% filter(sample != 'British_insample') 

# performance
tmp %>% 
  reshape2::dcast(trait + sample ~ method, value.var = 'r2_mean') %>%
  left_join(pop_map2, by = c('sample' = 'pop')) %>% mutate(pop2 = order_pop(pop2)) %>%
  ggplot() + 
  # geom_point(aes(x = PTRS, y = `(PT) PTRS`, color = '(PT) PTRS')) + 
  geom_point(aes(x = PTRS, y = `PTRS-b`, color = 'EN')) + 
  geom_point(aes(x = `(PT) PTRS`, y = `(PT) PTRS-b`, color = 'clump')) +
  facet_wrap(~pop2, scales = 'free') + 
  geom_abline(slope = 1, intercept = 0) +
  ylab('PTRS R2 with baboon genes') +
  xlab('PTRS R2 with all genes') + th2 +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# portability
tmp %>%
  left_join(pop_map2, by = c('sample' = 'pop')) %>% mutate(pop2 = order_pop(pop2)) %>%
  ggplot() +
  geom_violin(aes(x = pop2, y = portability, color = method), position = position_dodge(0.65)) +
  geom_boxplot(aes(x = pop2, y = portability, color = method), width = 0.1, position = position_dodge(0.65)) + th + 
  theme(legend.position = c(0.85, 0.85), legend.title = element_blank()) +
  theme(axis.title.x = element_blank())
