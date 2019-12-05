meta_info_to_query = function(trait_id, field_id, ninstance, narray) {
  out = list()
  for(i in 1 : length(trait_id)) {
    for(j in 1 : narray[i]) {
      for(n in 1 : ninstance[i]) {
        str = paste0('c', paste(field_id[i], n - 1, j - 1, sep = '_'))  # ninstance[i] - 1
        trait_name = paste0(trait_id[i], '-x-instance-', n - 1, '-x-array-', j - 1)  # ninstance[i] - 1
        out[[trait_name]] = str
      }
    }
  }
  out
}

trim_space = function(str) {
  stringr::str_remove_all(str, ' ')
}
