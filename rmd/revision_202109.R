library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size = 15))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')
score_color_code = c('PRS' = "#999999", 'PTRS (GTEx EUR)' = "#E69F00", 'PTRS (MESA EUR)' = "#E69F00", 'PTRS (MESA AFHI)' = "#56B4E9", 'PTRS (MESA ALL)' = "#F633FF")

df1 = load_perf('~/Desktop/tmp/ptrs-tf/from_nucleus/elastic_net_ptrs_gtex_british_updated_indivs_w_split_lambda.performance.csv') %>% filter(!trait %in% bad_traits)
df1pt = load_perf('~/Desktop/tmp/ptrs-tf/from_nucleus/pt_ptrs_gtex_british_updated_indivs_w_split_lambda.performance.csv') %>% filter(!trait %in% bad_traits)
df2 = load_perf('~/Desktop/tmp/ptrs-tf/from_nucleus/elastic_net_ptrs_mesa_british_updated_w_split_lambda.performance.csv') %>% filter(!trait %in% bad_traits)
df2_cau = df2 %>% filter(pred_expr_source == 'train')
df2_afhi = df2 %>% filter(pred_expr_source == 'against')
df3 = load_perf('~/Desktop/tmp/ptrs-tf/from_nucleus/elastic_net_ptrs_mesa_all_british_updated_w_split_lambda.performance.csv') %>% filter(!trait %in% bad_traits)
df4 = load_perf('~/Desktop/tmp/ptrs-tf/from_nucleus/prs_british_updated_indivs_w_split_lambda.performance.csv', is_prs = T) %>% filter(!trait %in% bad_traits)
df5 = load_perf('~/Desktop/tmp/ptrs-tf/from_nucleus/combine_british_updated_indivs_w_split_lambda.performance.csv', is_prs = T, max_n = 100000) %>% filter(!trait %in% bad_traits)
df5pt = load_perf('~/Desktop/tmp/ptrs-tf/from_nucleus/combine_pt_british_updated_indivs_w_split_lambda.performance.csv', is_prs = T, max_n = 100000) %>% filter(!trait %in% bad_traits)
df1s = get_test_perf_from_splits(df1)
df1spt = get_test_perf_from_splits(df1pt)
df2s_cau = get_test_perf_from_splits(df2_cau)
df2s_afhi = get_test_perf_from_splits(df2_afhi)
df3s = get_test_perf_from_splits(df3)
df4s = get_test_perf_from_splits(df4)
df5s = get_test_perf_from_splits(df5)
df5spt = get_test_perf_from_splits(df5pt)
df1ss = summarize_across_grps(df1s)
df1sspt = summarize_across_grps(df1spt)
df2ss_cau = summarize_across_grps(df2s_cau)
df2ss_afhi = summarize_across_grps(df2s_afhi)
df3ss = summarize_across_grps(df3s)
df4ss = summarize_across_grps(df4s)
df5ss = summarize_across_grps(df5s)
df5sspt = summarize_across_grps(df5spt)


df1ssp = calc_port(df1ss)
df1ssppt = calc_port(df1sspt)
df2ssp_cau = calc_port(df2ss_cau)
ref_mesa = df2ssp_cau %>% filter(sample == 'British_valid')
df2ssp_afhi = calc_port(df2ss_afhi, ref_mesa)
df3ssp = calc_port(df3ss)
df4ssp = calc_port(df4ss)
df5ssp = calc_port(df5ss)
df5ssppt = calc_port(df5sspt)

tmp = rbind(
  df1ssp %>% mutate(method = 'PTRS'),
  df1ssppt %>% mutate(method = '(PT) PTRS'),
  df5ssp %>% mutate(method = '(EN) PTRS+PRS'),
  df4ssp %>% mutate(method = 'PRS'),
  df5ssppt %>% mutate(method = '(PT) PTRS+PRS')
) %>% filter(sample != 'British_insample') 

tmp_mesa = rbind(
  df2ssp_cau %>% mutate(method = 'PTRS (MESA EUR)'),
  df2ssp_afhi %>% mutate(method = 'PTRS (MESA AFHI)'),
  df3ssp %>% mutate(method = 'PTRS (MESA ALL)'),
  df4ssp %>% mutate(method = 'PRS')
)

## generate files for figure 5 ###
{
  fig5_dir = 'analysis_output/combine_ptrs_prs/'
  dir.create(fig5_dir)
  p = tmp %>% 
    reshape2::dcast(trait + sample ~ method, value.var = 'r2_mean') %>% 
    ggplot() + 
    geom_point(aes(x = PRS, y = `(EN) PTRS+PRS`, color = 'PTRS+PRS')) + 
    geom_point(aes(x = PRS, y = `PTRS`, color = 'PTRS')) + 
    facet_wrap(~sample, scales = 'free') + 
    geom_abline(slope = 1, intercept = 0) +
    ylab('PTRS only/combine') + th2
  ggsave(paste0(fig5_dir, 'combine_ptrs_and_prs_en_r2.png'), p, width = 8, height = 5)
  
  p = tmp %>% 
    reshape2::dcast(trait + sample ~ method, value.var = 'r2_mean') %>% 
    ggplot() + 
    geom_point(aes(x = PRS, y = `(PT) PTRS+PRS`, color = 'PTRS+PRS')) + 
    geom_point(aes(x = PRS, y = `PTRS`, color = 'PTRS')) + 
    facet_wrap(~sample, scales = 'free') + 
    geom_abline(slope = 1, intercept = 0) +
    ylab('PTRS only/combine') + th2
  ggsave(paste0(fig5_dir, 'combine_ptrs_and_prs_pt_r2.png'), p, width = 8, height = 5)
  
  tmp %>% 
    reshape2::dcast(trait + sample ~ method, value.var = 'r2_mean') %>% 
    group_by(sample) %>% 
    do(test_a_vs_b(.$PRS, .$`(EN) PTRS+PRS`)) %>% 
    ungroup() %>%
    left_join(pop_map2, by = c('sample' = 'pop')) %>% 
    mutate(pop2 = order_pop(pop2)) %>% 
    select(pop2, delta, wilcox_pval) %>%
    rename(sample = pop2, `R2 (PRS - PTRS)` = delta) %>% 
    write.csv(paste0(fig5_dir, 'combine_ptrs_and_prs_en_r2.csv'))
  
  tmp %>% 
    reshape2::dcast(trait + sample ~ method, value.var = 'r2_mean') %>% 
    group_by(sample) %>% 
    do(test_a_vs_b(.$PRS, .$`(PT) PTRS+PRS`)) %>% 
    ungroup() %>%
    left_join(pop_map2, by = c('sample' = 'pop')) %>% 
    mutate(pop2 = order_pop(pop2)) %>% 
    select(pop2, delta, wilcox_pval) %>%
    rename(sample = pop2, `R2 (PRS - PTRS)` = delta) %>% 
    write.csv(paste0(fig5_dir, 'combine_ptrs_and_prs_pt_r2.csv'))
  
  p = tmp %>% filter(method %in% c('PRS', '(EN) PTRS+PRS')) %>% 
    # filter(sample == 'African') %>% 
    mutate(method = recode(method, `(EN) PTRS+PRS` = 'PTRS+PRS')) %>%
    left_join(pop_map2, by = c('sample' = 'pop')) %>% mutate(pop2 = order_pop(pop2)) %>%
    ggplot() +
    geom_violin(aes(x = pop2, y = portability, color = method), position = position_dodge(0.65)) +
    geom_boxplot(aes(x = pop2, y = portability, color = method), width = 0.1, position = position_dodge(0.65)) + th + 
    theme(legend.position = c(0.2, 0.15), legend.title = element_blank()) +
    theme(axis.title.x = element_blank())
  ggsave(paste0(fig5_dir, 'combine_ptrs_and_prs_en_portability.png'), p, width = 7, height = 5)
  
  p = tmp %>% filter(method %in% c('PRS', '(EN) PTRS+PRS')) %>% 
    filter(sample == 'African') %>%
    mutate(method = recode(method, `(EN) PTRS+PRS` = 'PTRS+PRS')) %>%
    left_join(pop_map2, by = c('sample' = 'pop')) %>% mutate(pop2 = order_pop(pop2)) %>%
    ggplot() +
    geom_violin(aes(x = method, y = portability, color = method), position = position_dodge(1)) +
    geom_boxplot(aes(x = method, y = portability, color = method), width = 0.2, position = position_dodge(1)) + th + 
    theme(legend.position = 'none', legend.title = element_blank()) +
    theme(axis.title.x = element_blank()); p
  ggsave(paste0(fig5_dir, 'combine_ptrs_and_prs_en_portability_zoom.png'), p, width = 4, height = 4)
  
  tmp %>% 
    reshape2::dcast(trait + sample ~ method, value.var = 'portability') %>% 
    group_by(sample) %>% 
    do(test_a_vs_b(.$PRS, .$`(EN) PTRS+PRS`)) %>% 
    ungroup() %>%
    left_join(pop_map2, by = c('sample' = 'pop')) %>% 
    mutate(pop2 = order_pop(pop2)) %>% 
    select(pop2, delta, wilcox_pval) %>%
    rename(sample = pop2, `portability (PRS - PTRS)` = delta) %>% 
    write.csv(paste0(fig5_dir, 'combine_ptrs_and_prs_en_portability.csv'))
  
  p = tmp %>% filter(method %in% c('PRS', '(PT) PTRS+PRS')) %>% 
    # filter(sample == 'African') %>% 
    mutate(method = recode(method, `(PT) PTRS+PRS` = 'PTRS+PRS')) %>%
    left_join(pop_map2, by = c('sample' = 'pop')) %>% mutate(pop2 = order_pop(pop2)) %>%
    ggplot() +
    geom_violin(aes(x = pop2, y = portability, color = method), position = position_dodge(0.65)) +
    geom_boxplot(aes(x = pop2, y = portability, color = method), width = 0.1, position = position_dodge(0.65)) + th + 
    theme(legend.position = c(0.2, 0.15), legend.title = element_blank()) +
    theme(axis.title.x = element_blank())
  ggsave(paste0(fig5_dir, 'combine_ptrs_and_prs_pt_portability.png'), p, width = 7, height = 5)
  
  p = tmp %>% filter(method %in% c('PRS', '(PT) PTRS+PRS')) %>% 
    filter(sample == 'African') %>%
    mutate(method = recode(method, `(PT) PTRS+PRS` = 'PTRS+PRS')) %>%
    left_join(pop_map2, by = c('sample' = 'pop')) %>% mutate(pop2 = order_pop(pop2)) %>%
    ggplot() +
    geom_violin(aes(x = method, y = portability, color = method), position = position_dodge(1)) +
    geom_boxplot(aes(x = method, y = portability, color = method), width = 0.2, position = position_dodge(1)) + th + 
    theme(legend.position = 'none', legend.title = element_blank()) +
    theme(axis.title.x = element_blank())
  ggsave(paste0(fig5_dir, 'combine_ptrs_and_prs_pt_portability_zoom.png'), p, width = 4, height = 4)
  
  tmp %>% 
    reshape2::dcast(trait + sample ~ method, value.var = 'portability') %>% 
    group_by(sample) %>% 
    do(test_a_vs_b(.$PRS, .$`(PT) PTRS+PRS`)) %>% 
    ungroup() %>%
    left_join(pop_map2, by = c('sample' = 'pop')) %>% 
    mutate(pop2 = order_pop(pop2)) %>% 
    select(pop2, delta, wilcox_pval) %>%
    rename(sample = pop2, `portability (PRS - PTRS)` = delta) %>% 
    write.csv(paste0(fig5_dir, 'combine_ptrs_and_prs_pt_portability.csv'))
}

## end ##


tmp %>% 
  reshape2::dcast(trait + sample ~ method, value.var = 'r2_mean') %>% 
  group_by(sample) %>% 
  do(test_a_vs_b(.$PRS, .$`(PT) PTRS+PRS`)) %>% pander::pander('R2 (PT) PTRS+PRS vs PRS')
tmp %>% 
  reshape2::dcast(trait + sample ~ method, value.var = 'r2_mean') %>% 
  group_by(sample) %>% 
  do(test_a_vs_b(.$PRS, .$`(EN) PTRS+PRS`)) %>% pander::pander('R2 (EN) PTRS+PRS vs PRS')
tmp %>% 
  reshape2::dcast(trait + sample ~ method, value.var = 'r2_mean') %>% 
  group_by(sample) %>% 
  do(test_a_vs_b(.$PRS, .$`PTRS`)) %>% pander::pander('R2 PTRS vs PRS')


tmp %>% 
  reshape2::dcast(trait + sample ~ method, value.var = 'portability') %>% 
  group_by(sample) %>% 
  do(test_a_vs_b(.$PRS, .$`(EN) PTRS+PRS`)) %>% 
  pander::pander('portability (EN) PTRS+PRS vs PRS')
tmp %>% 
  reshape2::dcast(trait + sample ~ method, value.var = 'portability') %>% 
  group_by(sample) %>% 
  do(test_a_vs_b(.$PRS, .$`PTRS`)) %>% 
  pander::pander('portability PTRS vs PRS')

tmp %>% filter(!method %in% c('PTRS', '(EN) PTRS+PRS')) %>% 
  filter(sample == 'African') %>% 
  left_join(pop_map2, by = c('sample' = 'pop')) %>% mutate(pop2 = order_pop(pop2)) %>%
  ggplot() +
  geom_violin(aes(x = pop2, y = portability, color = method), position = position_dodge(0.65)) +
  geom_boxplot(aes(x = pop2, y = portability, color = method), width = 0.1, position = position_dodge(0.65)) + th + 
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank())
tmp %>% 
  reshape2::dcast(trait + sample ~ method, value.var = 'portability') %>% 
  group_by(sample) %>% 
  do(test_a_vs_b(.$PRS, .$`(PT) PTRS`)) %>% 
  pander::pander('portability (PT) PTRS vs PRS')

tmp %>% filter(!method %in% c('PTRS', '(EN) PTRS+PRS')) %>% 
  # filter(sample == 'African') %>% 
  left_join(pop_map2, by = c('sample' = 'pop')) %>% mutate(pop2 = order_pop(pop2)) %>%
  ggplot() +
  geom_violin(aes(x = pop2, y = portability, color = method), position = position_dodge(0.65)) +
  geom_boxplot(aes(x = pop2, y = portability, color = method), width = 0.1, position = position_dodge(0.65)) + th + 
  theme(legend.position = c(0.15, 0.15), legend.title = element_blank())

outlier = tmp_mesa %>% filter(sample == 'African', portability > 0.5)
tmp_mesa %>% filter((!sample %in% 'African') | (!trait %in% outlier$trait)) %>% 
  filter(sample == 'African') %>% 
  left_join(pop_map2, by = c('sample' = 'pop')) %>% mutate(pop2 = order_pop(pop2)) %>%
  ggplot() +
  geom_violin(aes(x = pop2, y = portability, color = method), position = position_dodge(0.65)) +
  geom_boxplot(aes(x = pop2, y = portability, color = method), width = 0.1, position = position_dodge(0.65)) + th + 
  theme(legend.position = c(0.17, 0.8), legend.title = element_blank()) +
  scale_color_manual(values = score_color_code)

tmp_mesa %>% filter((!sample %in% 'African') | (!trait %in% outlier$trait)) %>% 
  filter(!sample %in% c('British_insample', 'British_valid')) %>% 
  left_join(pop_map2, by = c('sample' = 'pop')) %>% mutate(pop2 = order_pop(pop2)) %>%
  ggplot() +
  geom_violin(aes(x = pop2, y = portability, color = method), position = position_dodge(0.65)) +
  geom_boxplot(aes(x = pop2, y = portability, color = method), width = 0.1, position = position_dodge(0.65)) + th + 
  theme(legend.position = c(0.85, 0.8), legend.title = element_blank()) +
  scale_color_manual(values = score_color_code)

kk = rbind(
  tmp_mesa %>% filter((!sample %in% 'African') | (!trait %in% outlier$trait)) %>% 
    filter(sample == 'African') %>% 
    reshape2::dcast(trait + sample ~ method, value.var = 'portability') %>% 
    group_by(sample) %>% 
    do(test_a_vs_b(.$PRS, .$`PTRS (MESA EUR)`)) %>% 
    mutate(label = 'PRS vs PTRS (MESA EUR)'),
  tmp_mesa %>% filter((!sample %in% 'African') | (!trait %in% outlier$trait)) %>% 
    filter(sample == 'African') %>% 
    reshape2::dcast(trait + sample ~ method, value.var = 'portability') %>% 
    group_by(sample) %>% 
    do(test_a_vs_b(.$PRS, .$`PTRS (MESA AFHI)`)) %>% 
    mutate(label = 'PRS vs PTRS (MESA AFHI)'),
  tmp_mesa %>% filter((!sample %in% 'African') | (!trait %in% outlier$trait)) %>% 
    filter(sample == 'African') %>% 
    reshape2::dcast(trait + sample ~ method, value.var = 'portability') %>% 
    group_by(sample) %>% 
    do(test_a_vs_b(.$PRS, .$`PTRS (MESA ALL)`)) %>% 
    mutate(label = 'PRS vs PTRS (MESA ALL)')
)
