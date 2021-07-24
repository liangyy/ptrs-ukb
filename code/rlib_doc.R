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

load_hsq = function(fn) {
  df_h2 = read.delim2(fn)
  df_h2$h_sq = as.numeric(df_h2$h_sq)
  df_h2$h_sq_se = as.numeric(df_h2$h_sq_se)
  df_h2$h_sq[is.na(df_h2$h_sq)] = 0
  df_h2
}

order_pop = function(p) {
  factor(p, levels = c('EUR ref.', 'EUR test', 'EUR', 'S.ASN', 'E.ASN', 'CAR', 'AFR'))
}

load_perf = function(fn, ref_sample = 'British_valid', is_prs = F, simple = F) {
  df1 = read.csv(fn)
  if(isTRUE(is_prs)) {
    df1$sample[df1$sample == 'British_validation'] = 'British_valid'
    df1 = df1 %>% group_by(trait, sample, split_label) %>%
      mutate(model_name = paste0('model_', 0 : (n() - 1))) %>% ungroup()
  } else {
    df1 = df1 %>% group_by(trait, sample, split_label, pred_expr_source) %>% 
      mutate(model_name = paste0('model_', 0 : (n() - 1))) %>% ungroup()
  }
  if(isTRUE(simple)) {
    tmp = df1 %>% filter(sample == ref_sample) %>% filter(!is.na(partial_r2)) 
  } else {
    tmp = df1 %>% filter(split_label == 'repeat0_1')  %>% filter(sample == ref_sample) %>% filter(!is.na(partial_r2)) 
  }
  
  kk = list()
  for(tt in unique(tmp$trait)) {
    tmp2 = tmp %>% filter(trait == tt)
    if(nrow(tmp2) > 11) {
      tmp2 = tmp2[1 : 11, ]
    }
    kk[[length(kk) + 1]] = tmp2
  }
  tmp2 = do.call(rbind, kk)
  selected = paste0(tmp2$trait, tmp2$lambda)
  df1 = df1 %>% filter(paste0(trait, lambda) %in% selected)
}

get_test_perf_from_splits = function(df) {
  df = df %>% mutate(grp = unlist(lapply(strsplit(split_label, '_'), function(x) {x[1]})))
  get_split_r2 = function(clambda, csplit, cr2, tune_grp, eval_grp) {
    # partial_r2, split_label, lambda
    dd = data.frame(lambda = clambda, split_label = csplit, partial_r2 = cr2)
    tmp = dd %>% reshape2::dcast(lambda ~ split_label, value.var = 'partial_r2')
    x1 = tmp[[tune_grp]]
    x1[is.na(x1)] = -Inf
    x2 = tmp[[eval_grp]]
    x2[is.na(x2)] = -Inf
    list(value = x2[ which.max(x1)], model_name = tmp[['lambda']][ which.max(x1) ])
  }
  kk = df %>% group_by(trait, sample, alpha, grp) %>% 
    summarize(
      r2 = get_split_r2(model_name, split_label, partial_r2, tune_grp = paste0(grp[1], '_1'), eval_grp = paste0(grp[1], '_2'))$value, 
      model_name = get_split_r2(model_name, split_label, partial_r2, tune_grp = paste0(grp[1], '_1'), eval_grp = paste0(grp[1], '_2'))$model_name
    )
  kk
}

summarize_across_grps = function(dd) {
  dd_s = dd %>% group_by(trait, sample, alpha) %>% summarize(r2_mean = mean(r2), r2_sd = sd(r2), r2_median = median(r2), r2_qse = (quantile(r2, probs = 0.9) - quantile(r2, probs = 0.1)) / 2 / qnorm(0.9), model_name = paste0(unique(model_name), collapse = ',')) %>% ungroup() 
}

calc_port = function(dd, ref = NULL) {
  if(is.null(ref)) {
    ref = dd %>% filter(sample == 'British_valid')
  }
  others = dd # %>% filter(sample != 'British_valid')
  others %>% left_join(ref %>% select(trait, r2_mean, r2_sd), by = 'trait', suffix = c('', '.ref')) %>% mutate(portability = r2_mean / r2_mean.ref)
}

calc_port_max = function(dd, ref = NULL) {
  if(is.null(ref)) {
    ref = dd %>% filter(sample == 'British_valid')
  }
  others = dd # %>% filter(sample != 'British_valid')
  others %>% left_join(ref %>% select(trait, r2_max), by = 'trait', suffix = c('', '.ref')) %>% mutate(portability = r2_max / r2_max.ref)
}

calc_regu = function(dfh2, dfpve) {
  df_merge = inner_join(
    dfpve %>% select(trait, h_sq, h_sq_se),
    dfh2 %>% select(trait, h_sq, h_sq_se),
    by = 'trait', suffix = c('.pve', '.h2')
  )
  df_merge = df_merge %>% mutate(regu = h_sq.pve / h_sq.h2)
  ratio = delta_mtd(
    df_merge$h_sq.pve, df_merge$h_sq_se.pve ^ 2, 
    df_merge$h_sq.h2, df_merge$h_sq_se.h2 ^ 2
  )
  df_merge = df_merge %>% mutate(ratio_mean = ratio$m, ratio_se = sqrt(ratio$v))
  ratio_fe = meta_fixed(df_merge$ratio_mean, df_merge$ratio_se)
  list(fe = ratio_fe, raw = df_merge)
}

# test_afhi_vs_cau = function(scores, traits, labels) {
#   dd = data.frame(ss = scores, tt = traits, ll = labels) %>% 
#     reshape2::dcast(tt ~ ll, value.var = 'ss')
#   res = t.test(dd$AFHI, dd$CAU, paired = T)
#   diff = res$estimate
#   diff_95ci_low = res$conf.int[1]
#   diff_95ci_high = res$conf.int[2]
#   data.frame(delta_pve = diff, ci_low = diff_95ci_low, ci_high = diff_95ci_high)
# }
test_afhi_vs_cau = function(afhi, cau) {
  res = t.test(afhi, cau, paired = T)
  diff = res$estimate
  diff_95ci_low = res$conf.int[1]
  diff_95ci_high = res$conf.int[2]
  res2 = wilcox.test(afhi, cau, paired = TRUE)
  data.frame(delta_pve = diff, ci_low = diff_95ci_low, ci_high = diff_95ci_high, wilcox_stat = res2$statistic, wilcox_pval = res2$p.value)
}

test_a_vs_b = function(a, b) {
  res = t.test(a, b, paired = T)
  diff = res$estimate
  diff_95ci_low = res$conf.int[1]
  diff_95ci_high = res$conf.int[2]
  res2 = wilcox.test(a, b, paired = TRUE)
  data.frame(delta = diff, ci_low = diff_95ci_low, ci_high = diff_95ci_high, pval = res$p.value, wilcox_stat = res2$statistic, wilcox_pval = res2$p.value)
}

order_method = function(x) {
 factor(x, levels = c('PRS', 'PTRS (GTEx EUR)', 'PTRS (MESA EUR)', 'PTRS (MESA AFHI)', 'PTRS (MESA ALL)')) 
}

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
