meta_info_to_query = function(trait_id, field_id, ninstance, narray) {
  out = list()
  for(i in 1 : length(trait_id)) {
    for(j in 1 : narray[i]) {
      str = paste0('c', paste(field_id[i], ninstance[i] - 1, j - 1, sep = '_'))
      trait_name = paste0(trait_id[i], '-x-instance-', ninstance[i] - 1, '-x-array-', j - 1)
      out[[trait_name]] = str
    }
  }
  out
}
