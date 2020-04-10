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

trim_dot = function(str) {
  unlist(lapply(strsplit(str, '\\.'), function(x) { x[1] }))
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
  p = do.call(rbind, mydat) %>% ggplot() + geom_point(aes(x = x, y = y, color = color), ...) + facet_grid(cols = vars(xfacet), rows = vars(yfacet))
  p
}

parse_neale_snp = function(str) {
  unlist(lapply(strsplit(str, ':'), function(x) { paste0(x[1], ':', x[2]) }))
}

compute_r2 = function(df, y, ypred, covariates, report_pval = T) {
  covariates_terms = paste0(covariates, collapse = ' + ')
  formula_null = paste0('`', y, '`', ' ~ ', '1 + ', covariates_terms)
  formula_full = paste0('`', y, '`', ' ~ ', '1 + ', covariates_terms, ' + `', ypred, '`')
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

change_colname = function(df, from_name, to_name) {
  tmp = colnames(df)
  tmp[tmp == from_name] = to_name
  colnames(df) = tmp
  df
}

best_model_based_on_one = function(df, pop_name, model_col, score_col) {
  mydf = df %>% select(trait, population, subset)
  mydf$model = df[, model_col]
  mydf$score = df[, score_col]
  best_model = mydf %>% filter(population == pop_name) %>% group_by(trait, subset) %>% summarize(best_model = model[which.max(score)])
  perf_in_all = mydf %>% filter(paste(trait, subset, model) %in% paste(best_model$trait, best_model$subset, best_model$best_model)) %>% group_by(trait, subset) %>% mutate(transferability = score / score[population == pop_name])
  best_model = change_colname(best_model, 'model', model_col)
  best_model = change_colname(best_model, 'score', score_col)
  perf_in_all = change_colname(perf_in_all, 'model', model_col)
  perf_in_all = change_colname(perf_in_all, 'score', score_col)
  return(list(best_model = best_model, perf_in_all = perf_in_all))
}

best_model_for_each = function(df, reference_pop, model_col, score_col) {
  mydf = df %>% select(trait, population, subset)
  mydf$model = df[, model_col]
  mydf$score = df[, score_col]
  best_model = mydf %>% group_by(trait, subset, population) %>% summarize(best_model = model[which.max(score)])
  perf_in_all = mydf %>% filter(paste(trait, subset, model, population) %in% paste(best_model$trait, best_model$subset, best_model$best_model, best_model$population)) %>% group_by(trait, subset) %>% mutate(transferability = score / score[population == reference_pop])
  best_model = change_colname(best_model, 'model', model_col)
  best_model = change_colname(best_model, 'score', score_col)
  perf_in_all = change_colname(perf_in_all, 'model', model_col)
  perf_in_all = change_colname(perf_in_all, 'score', score_col)
  return(list(best_model = best_model, perf_in_all = perf_in_all))
}

delta_mtd = function(mx, vx, my, vy) {
  m = mx / my
  v = vx / (my ^ 2) + vy * (mx ^ 2) / (my ^ 4)
  return(list(m = m, v = v))
}
meta_fixed = function(m, se) {
  keep_ind = !is.na(m) & !is.infinite(m) & !is.nan(m) & se != 0 & !is.na(se)
  m = m[keep_ind]
  se = se[keep_ind]
  w = 1 / (se^2)
  mfe = sum(m * w) / sum(w)
  vfe = 1 / sum(w)
  return(list(m = mfe, se = sqrt(vfe)))
}

get_meta_for_supp = function() {
  color_category = list()
  color_category[['Blood cell counts']] = c('wbc', 'rbc', 'platelet', 'lymphocyte', 'monocyte', 'neutrophil', 'eosinophil', 'basophil')
  color_category[['Haemoglobin related']] = c('mcv', 'mch', 'mchc', 'hb', 'ht')
  color_category[['Blood pressures']] = c('dbp', 'sbp')
  color_category[['Height']] = c('height')
  color_category[['BMI']] = c('bmi')
  color_mixer_candidates = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # c('blood_cell_counts' = '#009E73', 'haemoglobin_related' = '#F0E442', 'blood_pressures' = '#0072B2', 'height' = '#D55E00', 'bmi' = '')
  color_mixer = color_mixer_candidates[c(-1, -2, -3)]
  names(color_mixer) = names(color_category)
  color_category_guide = c()
  color_category_guide_name = c()
  for(n in names(color_category)) {
    for(l in color_category[[n]]) {
      color_category_guide_name = c(color_category_guide_name, l)
      color_category_guide = c(color_category_guide, n)
    }
  }
  df_color_category = data.frame(trait = color_category_guide_name, group = color_category_guide)
  pop_color_mixer = c('African' = '#E74C3C', 'British_test' = '#28B463', 'British_validation' = '#D4AC0D', 'Chinese' = '#3498DB', 'Indian' = '#9B59B6', 'British_insample' = '#AAB7B8', 'British \n (test)' = '#28B463', 'British \n (validation)' = '#D4AC0D', 'British \n (in sample)' = '#AAB7B8')
  type_shape = c('PTRS' = 0, 'PRS' = 4)
  score_color_code = c("PRS" = "#999999", "PTRS (GTEx EUR)" = "#E69F00", "PTRS (MESA CAU)" = "#E69F00", "PTRS (MESA CAU & AFHI)" = "#56B4E9")
  return(list(df_color_category = df_color_category, color_mixer = color_mixer, pop_color_mixer = pop_color_mixer, type_shape = type_shape, score_color_code = score_color_code))
}

update_popname = function(p, type = 1) {
  if(type == 1) {
    p[p == 'British_insample'] = 'British \n (in sample)'
    p[p == 'British_test'] = 'British \n (test)'
    p[p == 'British_validation'] = 'British \n (validation)'
    p = factor(p, levels = c('British \n (in sample)', 'British \n (validation)', 'British \n (test)', 'Indian', 'Chinese', 'African'))
  } else if(type == 2) {
    p[p == 'British_insample'] = 'British'
    p[p == 'British_test'] = 'British'
    p[p == 'British_validation'] = 'British'
    p = factor(p, levels = c('British', 'Indian', 'Chinese', 'African'))
  }
  p
}

order_score = function(p) {
  factor(p, levels = c('PRS', 'PTRS (MESA CAU)', 'PTRS (MESA CAU & AFHI)'))
}

# 
# my_parser = function(x) {
#   as.numeric(stringr::str_remove(x, 'pval_'))
# }
