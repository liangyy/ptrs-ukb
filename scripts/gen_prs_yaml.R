nsubset = 17
traits = tolower(read.csv('../external_data/martin_et_al_2019ng_table_s6.csv')$Trait)
other_populations = c('African', 'Chinese', 'Indian')
indiv_list_prefix = '/vol/bmd/yanyul/GitHub/ptrs-ukb/output/data_split/'
sum_stat_prefix = '/vol/bmd/yanyul/UKB/gwas_on_subset/gwas_runs_in_tsv/gwas_in_tsv_'
ld_clump_prefix = '/vol/bmd/yanyul/UKB/ld_clump/gwas_clump_x_'
out_list = list()
for(i in 1 : nsubset) {
  list_name = paste0('subset', i)
  out_list[[list_name]] = list()
  indiv_lists = c(
    other_populations, 
    c(
      paste0('British-test-', i, '.txt'),
      paste0('British-validation-', i, '.txt')
    )
  )
  indiv_lists_path = paste0(indiv_list_prefix, indiv_lists)
  out_list[[list_name]][['indiv_lists']] = paste0(indiv_lists_path, collapse = ',')
  out_list[[list_name]][['GWASs']] = list()
  for(t in traits) {
    out_list[[list_name]][['GWASs']][[t]] = list()
    out_list[[list_name]][['GWASs']][[t]][['sum_stat']] = paste0(sum_stat_prefix, list_name, '_x_', t, '.tsv')
    out_list[[list_name]][['GWASs']][[t]][['ld_clump']] = paste0(ld_clump_prefix, list_name, '_x_', t, '.clumped_snp')
  }
}

yaml::write_yaml(out_list, '../misc/prs_full.yaml')
