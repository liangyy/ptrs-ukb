meta_info_to_query = function(trait_id, field_id, ninstance, narray) {
  out = list()
  for(i in 1 : length(trait_id)) {
    for(j in 1 : narray[i]) {
      for(n in 1 : ninstance[i]) {
        str = paste0('c', paste(field_id[i], n - 1, j - 1, sep = '_'))  # ninstance[i] - 1
        trait_name = paste0(trait_id[i], '_x_instance_', n - 1, '_x_array_', j - 1)  # ninstance[i] - 1
        out[[trait_name]] = str
      }
    }
  }
  out
}

trim_space = function(str) {
  stringr::str_remove_all(str, ' ')
}

first_non_na = function(vec) {
  out = NA
  if(sum(!is.na(vec)) > 0) {
    out = vec[which(!is.na(vec))[1]]
  }
  out
}

take_average = function(vec) {
  out = NA
  if(sum(!is.na(vec)) > 0) {
    out = mean(vec, na.rm = T)
  }
  out
}

aggregate_instances = function(df, method = first_non_na) {
  apply(df, 1, first_non_na)
}

sample_at_most_n = function(n, nmax) {
  if(n > nmax) {
    o = 1 : n
    o = o %in% sample(1 : n, size = nmax, replace = F)
    return(o)
  } else {
    return(rep(TRUE, n))
  }
}

subsample_for_vis = function(group, nmax = 10000) {
  e = data.frame(group = group)
  e = e %>% group_by(group) %>% mutate(is_selected = sample_at_most_n(n(), nmax = nmax))
  e %>% pull(is_selected)
}

myggpairs = function(df, col, ...) {
  mydat = list()
  done = c()
  for(i in colnames(df)) {
    for(j in colnames(df)) {
      if(i == j) {
        next
      }
      if(paste(i, j) %in% done) {
        next
      }
      tmp = df[, c(i, j)]
      colnames(tmp) = c('x', 'y')
      mydat[[length(mydat) + 1]] = tmp %>% mutate(color = col, xfacet = i, yfacet = j)
      done = c(done, paste(i, j), paste(j, i))
    }
  }
  p = do.call(rbind, mydat) %>% ggplot() + geom_point(aes(x = x, y = y, color = color), ...) + facet_grid(xfacet ~ yfacet)
  p
}

parse_neale_snp = function(str) {
  unlist(lapply(strsplit(str, ':'), function(x) { paste0(x[1], ':', x[2]) }))
}

compute_r2 = function(df, y, ypred, covariates, report_pval = T) {
  covariates_terms = paste0(covariates, collapse = ' + ')
  formula_null = paste0(y, ' ~ ', '1 + ', covariates_terms)
  formula_full = paste0(y, ' ~ ', '1 + ', covariates_terms, ' + `', ypred, '`')
  # print(formula_null)
  # print(formula_full)
  mod_null <- lm(as.formula(formula_null), data = df)
  mod_full <- lm(as.formula(formula_full), data = df)
  # derived from asbio::partial.R2
  SSE.wo <- get_sse(mod_null, df[, y])
  SSE.with <- get_sse(mod_full, df[, y])
  r2 <- (SSE.wo - SSE.with) / SSE.wo
  if(report_pval == F) {
    return(data.frame(r2 = r2, SSE.wo = SSE.wo, SSE.with = SSE.with))
  } else {
    pval <- anova(mod_full, mod_null)$'Pr(>F'[2]
    return(data.frame(r2 = r2, pval = pval, SSE.wo = SSE.wo, SSE.with = SSE.with))
  }
}

get_sse = function(mod, y) {
  ypred = predict(mod)
  sum((y - ypred)^2)
}

report_r2 = function(df, y, ypred, covariates, nbootstrap = 1000, quantiles = c(0.025, 0.975)) {
  obs = compute_r2(df, y, ypred, covariates)
  bootstrap_out = replicate(nbootstrap, {
    shuffled_idx = sample(1 : nrow(df), size = nrow(df), replace = T)
    compute_r2(df[shuffled_idx, ], y, ypred, covariates, report_pval = F)$r2
  })
  r2_quantile = quantile(bootstrap_out, probs = quantiles)
  r2_quantile = as.data.frame(t(r2_quantile))
  colnames(r2_quantile) = paste0('quantile_', quantiles)
  rownames(r2_quantile) = NULL
  cbind(obs, r2_quantile)
}
