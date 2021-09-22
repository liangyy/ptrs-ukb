library(ggplot2)
library(dplyr)
theme_set(theme_classic(base_size = 15))

pop_map = data.frame(
  pop = c('African', 'Chinese', 'Indian', 'Caribbean', 'British_test', 'British_valid'),
  pop2 = c('AFR', 'E.ASN', 'S.ASN', 'CAR', 'EUR (test)', 'EUR (validation)')
)

# comparing previous version of R2 vs the new definition
load_best = function(fn, desired_tag = NULL, mode = 'ptrs') {
  if(mode == 'ptrs') {
    df_en = read.csv(fn)
    if(!is.null(desired_tag)) {
      df_en = df_en %>% filter(pred_expr_source == desired_tag)
    }
    # add model name
    df_en = df_en %>% group_by(trait, alpha, sample) %>% 
      mutate(model_name = paste0('model_', 0 : (n() - 1))) %>% 
      ungroup()
    df_en = df_en %>% 
      filter(!is.na(partial_r2)) %>% 
      group_by(trait, alpha, sample) %>% 
      mutate(order = 1 : n()) %>% 
      ungroup() %>% 
      filter(order <= 11) %>% 
      select(-order)
    best = df_en %>% select(sample, trait, partial_r2, model_name) %>%
      group_by(trait, sample) %>% 
      summarize(
        r2_max = max(partial_r2, na.rm = T), 
        model_name = model_name[which.max(partial_r2)]
      ) %>% ungroup()
  } else if(mode == 'prs') {
    df_prs = read.table(fn, header = T, sep = '\t', stringsAsFactors = F)
    df_prs$sample[df_prs$sample == 'British-test-1'] = 'British_test'
    df_prs$sample[df_prs$sample == 'British-validation-1'] = 'British_valid'
    df_prs = df_prs %>% rename(prs_cutoff = ptrs_cutoff)
    df_prs = df_prs %>% arrange(desc(-prs_cutoff))
    # add model name
    df_prs = df_prs %>% group_by(trait, sample) %>% 
      mutate(model_name = paste0('model_', 0 : (n() - 1))) %>% 
      ungroup()
    best = df_prs %>% group_by(trait, sample) %>% 
      summarize(
        r2_max = max(partial_r2, na.rm = T), 
        model_name = model_name[which.max(partial_r2)]
      ) %>% ungroup()
  }
  best
}
load_r2 = function() {
  df_prev = rbind(
    load_best('../ptrs-analysis/data/elastic_net_ptrs_gtex_british.performance.csv') %>% select(r2_max, trait, sample, model_name) %>% mutate(tag = 'GTEx') ,
    load_best('../ptrs-analysis/data/elastic_net_ptrs_mesa_british.performance.csv', desired_tag = 'against') %>% select(r2_max, trait, sample, model_name)  %>% mutate(tag = 'MESA_AFHI'),
    load_best('../ptrs-analysis/data/elastic_net_ptrs_mesa_british.performance.csv', desired_tag = 'train') %>% select(r2_max, trait, sample, model_name)  %>% mutate(tag = 'MESA_EUR'),
    load_best('../ptrs-analysis/data/partial_r2-prs.subset1_British.tsv', mode = 'prs') %>% select(r2_max, trait, sample, model_name)  %>% mutate(tag = 'PRS')
  ) %>% rename(R2 = r2_max)
  df_new = readRDS('analysis_output/revision_june_2021.r2.rds')
  tmp = list()
  for(i in names(df_new)) {
    tmp[[i]] = df_new[[i]] %>% select(r2_mean, trait, sample, model_name) %>% mutate(tag = i) %>% rename(R2 = r2_mean)
  }
  df_new = do.call(rbind, tmp)
  return(list(prev = df_prev, new = df_new))
}
load_genes = function(dn, traits) {
  get_idx = function(ss) {
    as.numeric(strsplit(ss, '_')[[1]][2])
  }
  res = list()
  for(trait in traits) {
    res[[trait]] = list()
    fn = paste0(dn, '/weights.', trait, '.tsv.gz')
    tmp = read.table(fn, header = T)
    # if all are not degenerative, the max idx is 19
    offset = 19 - get_idx(colnames(tmp)[ncol(tmp)])
    for(col in colnames(tmp)[-1]) {
      idx = get_idx(col)
      new_col = paste0('model_', idx + offset)
      res[[trait]][[new_col]] = tmp$gene_id[abs(tmp[[col]]) > 1e-3]
    }
  }
  res
}
gen_upsetr = function(df, gene_sets) {
  trait_list = list()
  for(tt in unique(df$trait)) {
    tmp = df %>% filter(trait == tt)
    gene_list_inner = list()
    gene_all = c()
    for(ss in unique(tmp$sample)) {
      tmp2 = tmp %>% filter(sample == ss)
      if(nrow(tmp2) != 1) {
        message('Something wrong.')
        next
      }
      models = strsplit(tmp2$model_name[1], ',')[[1]]
      genes = c()
      for(model in models) {
        genes = union(genes, unique(gene_sets[[tt]][[model]]))
      }
      gene_all = union(gene_all, genes)
      gene_list_inner[[ss]] = genes
    }
    gene_mat = data.frame(gene = gene_all)
    for(ss in names(gene_list_inner)) {
      curr = gene_mat$gene %in% gene_list_inner[[ss]]
      gene_mat[[ss]] = as.numeric(curr)
    }
    trait_list[[tt]] = gene_mat
  }
  trait_list
}
plot_upset = function(res, dn) {
  dir.create(paste0('analysis_output/', dn))
  for(tt in names(res)) {
    png(file = paste0('analysis_output/', dn, '/', tt, '.png'), 
        res = 300, width = 4.75, height = 3.25, units="in")
    cols = colnames(res[[tt]])[-1]
    cols = data.frame(pop = cols) %>% left_join(pop_map, by = 'pop')
    tmp = res[[tt]]
    colnames(tmp)[-1] = cols$pop2
    tmp = tmp[, !(is.na(colnames(tmp)) | colnames(tmp) %in% c('EUR (validation)'))]
    print(UpSetR::upset(tmp))
    print(grid::grid.text(tt,x = 0.65, y=0.95, gp=grid::gpar(fontsize = 8)))
    dev.off()
  }
}

# main
res = load_r2()
df_prev = res$prev
df_new = res$new
df_merge = inner_join(df_prev, df_new, by = c('trait', 'sample', 'tag'), suffix = c('.prev', '.new'))
table(df_merge$sample, df_merge$tag)
# filter out British insample and Indian
# since PRS does not have British insample result
# and the sample size of Indian is changed so the 
# results are not comparable
df_merge = df_merge %>% filter(! sample %in% c('British_insample', 'Indian'))
pp = df_merge %>% 
  left_join(pop_map, by = c('sample' = 'pop')) %>%
  ggplot() +
  geom_abline(slope = 1, intercept = 0) + 
  geom_point(aes(x = R2.prev, y = R2.new, color = tag), alpha = 0.5, size = 5) + 
  facet_wrap(~pop2, scales = 'free') 
ggsave('analysis_output/revision_misc.r2.png', pp)


genes = load_genes('~/Desktop/tmp/ptrs-tf/from_nucleus/models_new/elastic_net_ptrs_gtex_british.revision_0.1.export_model', unique(df_merge$trait))
res = gen_upsetr(df_new %>% filter(tag == 'GTEx'), genes)
dn = 'revision_misc_upsetr_gtex'
plot_upset(res, dn)

genes = load_genes('~/Desktop/tmp/ptrs-tf/from_nucleus/models_new/elastic_net_ptrs_mesa_cau_british.revision_0.1.export_model', unique(df_merge$trait))
res = gen_upsetr(df_new %>% filter(tag == 'MESA_EUR'), genes)
dn = 'revision_misc_upsetr_mesa_eur'
plot_upset(res, dn)

genes = load_genes('~/Desktop/tmp/ptrs-tf/from_nucleus/models_new/elastic_net_ptrs_mesa_all_british.revision_0.1.export_model', unique(df_merge$trait))
res = gen_upsetr(df_new %>% filter(tag == 'MESA_ALL'), genes)
dn = 'revision_misc_upsetr_mesa_all'
plot_upset(res, dn)

# portability using no splitting
df_ptrs = load_best('~/Desktop/tmp/ptrs-tf/from_nucleus/elastic_net_ptrs_gtex_british_updated_indivs_lambda.performance.csv') %>% select(r2_max, trait, sample, model_name) %>% mutate(tag = 'GTEx')
df_prs = load_best('~/Desktop/tmp/ptrs-tf/from_nucleus/prs_british_updated_indivs_lambda.performance.csv') %>% select(r2_max, trait, sample, model_name) %>% mutate(tag = 'PRS')
df_prs$sample[ df_prs$sample == 'British_validation' ] = 'British_valid'
df_merge = rbind(df_prs, df_ptrs %>% filter(sample != 'British_insample'))
df_merge = df_merge %>% group_by(trait, tag) %>% mutate(portability = r2_max / r2_max[sample == 'British_valid']) %>% ungroup()
df_merge %>% ggplot() + 
  geom_violin(aes(x = sample, y = portability, color = tag), position = position_dodge(0.65)) +
  geom_boxplot(aes(x = sample, y = portability, color = tag), width = 0.1, position = position_dodge(0.65)) + 
  th + 
  # scale_color_manual(values = score_color_code) + 
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank())
df_merge %>% 
  reshape2::dcast(trait + sample ~ tag, value.var = 'portability') %>% 
  group_by(sample) %>% 
  do(test_a_vs_b(.$PRS, .$GTEx)) %>% 
  pander('Portability: GTEx EUR whole blood')

df_ptrs_mesa_cau = load_best('~/Desktop/tmp/ptrs-tf/from_nucleus/elastic_net_ptrs_mesa_british_updated_lambda.performance.csv', desired_tag = 'train') %>% select(r2_max, trait, sample, model_name) %>% mutate(tag = 'MESA_EUR')
df_ptrs_mesa_afhi = load_best('~/Desktop/tmp/ptrs-tf/from_nucleus/elastic_net_ptrs_mesa_british_updated_lambda.performance.csv', desired_tag = 'against') %>% select(r2_max, trait, sample, model_name) %>% mutate(tag = 'MESA_AFHI')
df_merge_mesa_cau = rbind(df_prs, df_ptrs_mesa_cau %>% filter(sample != 'British_insample')) %>% 
  group_by(trait, tag) %>% 
  mutate(portability = r2_max / r2_max[sample == 'British_valid']) %>% 
  ungroup()
tmp = rbind(
  df_prs, df_ptrs_mesa_cau, df_ptrs_mesa_afhi) %>% 
  filter(sample != 'British_insample')
df_merge_mesa_cau %>% ggplot() + 
  geom_violin(aes(x = sample, y = portability, color = tag), position = position_dodge(0.65)) +
  geom_boxplot(aes(x = sample, y = portability, color = tag), width = 0.1, position = position_dodge(0.65)) + 
  th + 
  # scale_color_manual(values = score_color_code) + 
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank())
kk_cau = calc_port(df_ptrs_mesa_cau %>% rename(r2_mean = r2_max) %>% mutate(r2_sd = 0))
kk_prs = calc_port(df_prs %>% rename(r2_mean = r2_max) %>% mutate(r2_sd = 0))
kk_afhi = calc_port(
  df_ptrs_mesa_afhi %>% rename(r2_mean = r2_max) %>% mutate(r2_sd = 0),
  df_ptrs_mesa_cau %>% filter(sample == 'British_valid') %>% rename(r2_mean = r2_max) %>% mutate(r2_sd = 0))
rbind(
  kk_afhi, kk_cau, kk_prs) %>% filter(sample == 'African') %>% 
  ggplot() + 
  geom_violin(aes(x = sample, y = portability, color = tag), position = position_dodge(0.65)) +
  geom_boxplot(aes(x = sample, y = portability, color = tag), width = 0.1, position = position_dodge(0.65)) + 
  th

