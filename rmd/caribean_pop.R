# dd = data.table::fread('output/new_query_phenotypes_cleaned_up.csv', sep = ',', data.table = F)
# old_list = rbind(
#   read.table('output/data_split/African.txt', header = T) %>% mutate(pop = 'African'),
#   read.table('output/data_split/Indian.txt', header = T) %>% mutate(pop = 'Indian')
# )
# new_list = rbind(
#   read.table('output/new_data_split/African.txt', header = T) %>% mutate(pop = 'African'),
#   read.table('output/new_data_split/Indian.txt', header = T) %>% mutate(pop = 'Indian')
# )
# eur_list = read.table('output/data_split/British-test-1.txt', header = T) %>% mutate(pop = 'British')
# ch_list = read.table('output/data_split/Chinese.txt', header = T) %>% mutate(pop = 'Chinese')
# 
# dd = dd[ dd$eid %in% c(old_list$IID, new_list$IID, eur_list$IID, ch_list$IID), ]
# indiv_label = rep(NA, nrow(dd))
# indiv_label[ dd$eid %in% eur_list$IID ] = 'EUR'
# indiv_label[ dd$eid %in% ch_list$IID ] = 'E.ASN'
# indiv_label[ dd$eid %in% old_list$IID[old_list$pop == 'Indian'] ] = 'old-S.ASN'
# indiv_label[ dd$eid %in% old_list$IID[old_list$pop == 'African'] ] = 'old-AFR'
# indiv_label[is.na(indiv_label) & dd$eid %in% new_list$IID[new_list$pop == 'Indian'] ] = 'new-S.ASN'
# indiv_label[is.na(indiv_label) & dd$eid %in% new_list$IID[new_list$pop == 'African'] ] = 'new-AFR'
# sum(is.na(indiv_label))
# dd$pop_label = indiv_label
# dd %>% ggplot() + geom_point(aes(x = pc1, y = pc2, color = indiv_label))
# dd %>% ggplot() + geom_point(aes(x = pc2, y = pc3, color = indiv_label))
# dd %>% ggplot() + geom_point(aes(x = pc1, y = pc3, color = indiv_label))
# dd %>% ggplot() + geom_point(aes(x = pc1, y = pc4, color = indiv_label))
# dd %>% ggplot() + geom_point(aes(x = pc2, y = pc4, color = indiv_label))
# dd %>% ggplot() + geom_point(aes(x = pc3, y = pc4, color = indiv_label))
# dd %>% ggplot() + geom_point(aes(x = pc3, y = pc4, color = indiv_label))
# 
# kk = data.table::fread('output/new_query_phenotypes_output.csv', sep = ',', data.table = F)
# head(kk)
# mm = kk[, c('eid', 'ethnicity_x_instance_0_x_array_0', 'ethnicity_x_instance_1_x_array_0', 'ethnicity_x_instance_2_x_array_0')]
# mm_multi = mm[rowSums(is.na(mm)) > 3, ]

df = data.table::fread('../output/new2_query_phenotypes_cleaned_up.csv', sep = ',', data.table = F)
df0 = data.table::fread('../output/new_query_phenotypes_cleaned_up.csv', sep = ',', data.table = F)
df0 = left_join(df0, df %>% select(eid, meaning), by = 'eid', suffix = c('.old', ''))
df0 = left_join(df0, dd$dat_ll %>% select(eid, `African dist`), by = 'eid')
eur_list = read.table('../output/data_split/British-test-1.txt', header = T) %>% mutate(pop = 'British')
d1 = df0 %>% filter((meaning == 'British' & eid %in% eur_list$IID) | meaning != 'British') 
# # df %>% filter(meaning %in% c('African', 'Caribbean'))
# d1 %>% ggplot() + geom_boxplot(aes(x = meaning, y = height)) + facet_wrap(~sex)
# d1 %>% ggplot() + geom_point(aes(x = pc4, y = height, color = meaning), alpha = 0.1) + facet_grid(meaning~sex)
# d1 %>% ggplot() + geom_histogram(aes(x = pc4, color = meaning), position = 'dodge')
d1 %>% ggplot() + geom_point(aes(x = pc1, y = pc2, color = meaning), alpha = 0.3)
d1 %>% ggplot() + geom_point(aes(x = pc2, y = pc3, color = meaning), alpha = 0.3)
d1 %>% ggplot() + geom_point(aes(x = pc1, y = pc3, color = meaning), alpha = 0.3)
d1 %>% ggplot() + geom_point(aes(x = pc1, y = pc4, color = meaning), alpha = 0.3)
d1 %>% ggplot() + geom_point(aes(x = pc2, y = pc4, color = meaning), alpha = 0.3)
d1 %>% ggplot() + geom_point(aes(x = pc3, y = pc4, color = meaning), alpha = 0.3)
# d1 %>% ggplot() + geom_point(aes(x = pc3, y = pc4, color = meaning), alpha = 0.3)
# 
# d1 %>% ggplot() + geom_point(aes(x = pc3, y = pc4, color = meaning))
# cores = d1 %>% filter(pc3 <= 0, pc3 >= -12.5, pc4 >= -12.5, pc4 <= 12.5)
# cores %>% ggplot() + geom_point(aes(x = pc1, y = pc2, color = meaning), alpha = 0.1)
# cores %>% ggplot() + geom_point(aes(x = pc2, y = pc3, color = meaning), alpha = 0.1)
# cores %>% ggplot() + geom_point(aes(x = pc1, y = pc3, color = meaning), alpha = 0.1)
# cores %>% ggplot() + geom_point(aes(x = pc1, y = pc4, color = meaning), alpha = 0.1)
# cores %>% ggplot() + geom_point(aes(x = pc2, y = pc4, color = meaning), alpha = 0.1)
# cores %>% ggplot() + geom_point(aes(x = pc3, y = pc4, color = meaning), alpha = 0.1)
# cores2 = cores %>% filter(pc1 > 350)
# 
# d1[, c('meaning', paste0('pc', 1:20)) ] %>% reshape2::melt() %>% ggplot() + geom_histogram(aes(x = value, color = meaning), position = 'dodge') + facet_wrap(~variable, scales = 'free')
