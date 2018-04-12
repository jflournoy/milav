#' check_manifests_exist
#'
#' @param manifest_vec
#' @param aDF
#'
#' @return
#' @export
check_manifests_exist <- function(manifest_vec, aDF){
  df_names <- names(aDF)
  wave_manifest_names <- unlist(lapply(c('a', 'b', 'c', 'd'), paste0, manifest_vec))
  all_in_df_names <- all(wave_manifest_names %in% df_names)
  return(all_in_df_names)
}

#' create_lav_id_strat_widaman
#'
#' @param factor_name
#'
#' @return
#' @export
create_lav_id_strat_widaman <- function(factor_name){
  factor_name_w1 <- paste0(factor_name, '_W1')
  factor_intercept <- paste0(factor_name_w1, ' ~ 0*1')
  factor_var <- paste0(factor_name_w1, ' ~~ 1*', factor_name_w1)
  paste0(factor_intercept, '\n', factor_var)
}

#' create_lav_growth_model
#'
#' @param factor_name
#' @param manifest_vec
#' @param int.i.group_equal
#' @param var.i.group_equal
#' @param int.s.group_equal
#' @param var.s.group_equal
#' @param cov.group_equal
#' @param num_groups
#' @param center_wave
#'
#' @return
#' @export
create_lav_growth_model <- function(factor_name, manifest_vec, int.i.group_equal = F, var.i.group_equal = F, int.s.group_equal = F, var.s.group_equal = F, cov.group_equal = F, num_groups = 4, center_wave = 1){
  wave_factors <- paste0(factor_name, '_W', 1:4)
  iname <- paste0(factor_name, '_i')
  intercept <- paste0(iname, ' =~ ',
                      paste(
                        paste0('1*', wave_factors),
                        collapse = ' + '))
  sname <- paste0(factor_name, '_s')
  slope <- paste0(sname, ' =~ ',
                  paste(
                    paste0(1:4 - center_wave, '*', wave_factors),
                    collapse = ' + '))
  factor_intercepts <- paste(paste0(wave_factors, ' ~ 0*1'), collapse = '\n')
  factor_variances <- paste(
    paste0(wave_factors,
           paste0('~~ c(', rep('fctrvar', num_groups), ')*'),
           wave_factors),
    collapse = '\n')
  growth_code <- paste(intercept, slope, factor_intercepts, factor_variances, sep = '\n')
  if(int.i.group_equal){
    int.i <- paste0(iname, ' ~ c(', paste(rep('inti', num_groups), collapse = ', '), ')*1')
    growth_code <- paste(growth_code, int.i, sep = '\n')
  }
  if(int.s.group_equal){
    int.s <- paste0(sname, ' ~ c(', paste(rep('ints', num_groups), collapse = ', '), ')*1')
    growth_code <- paste(growth_code, int.s, sep = '\n')
  }
  if(var.i.group_equal){
    var.i <- paste0(iname, ' ~~ c(', paste(rep('vari', num_groups), collapse = ', '), ')*', iname)
    growth_code <- paste(growth_code, var.i, sep = '\n')
  }
  if(var.s.group_equal){
    var.s <- paste0(sname, ' ~~ c(', paste(rep('vars', num_groups), collapse = ', '), ')*', sname)
    growth_code <- paste(growth_code, var.s, sep = '\n')
  }
  if(cov.group_equal){
    covar.is <- paste0(sname, ' ~~ c(', paste(rep('covar_is', num_groups), collapse = ', '), ')*', iname)
    growth_code <- paste(growth_code, covar.is, sep = '\n')
  }
  return(growth_code)
}

#' create_lav_factor_loadings
#'
#' @param factor_name
#' @param manifest_vec
#' @param group_equal
#' @param long_free
#' @param num_groups
#'
#' @return
#' @export
create_lav_factor_loadings <- function(factor_name, manifest_vec, group_equal = F, long_free = F, num_groups = 4){
  if(long_free && group_equal){
    stop("Can't set groups equal and free longitudinal constraints.")
  }
  wave_factors <- paste0(factor_name, '_W', 1:4)
  num_indicators <- length(manifest_vec)
  load_labels <- lapply(paste0('L', 1:num_indicators),
                        function(label){
                          rep_labs <- rep(label, num_groups)
                          if(!group_equal){
                            rep_labs <- paste0(rep_labs, letters[1:num_groups])
                          }
                          paste0('c(', paste(rep_labs, collapse = ', '), ')')
                        })
  loadings <- c(paste0(load_labels, '*'))
  wave_manifests_rhs <- lapply(letters[1:4],
                               function(letter){
                                 if(long_free) {
                                   paste(paste0(c(loadings[[1]], rep('', length(loadings)-1)),
                                                letter,
                                                manifest_vec),
                                         collapse = ' + ')
                                 } else {
                                   paste(paste0(loadings,
                                                letter,
                                                manifest_vec),
                                         collapse = ' + ')
                                 }
                               })
  wave_eqs <- paste0(wave_factors, ' =~ ', wave_manifests_rhs)
  all_wave_model <- paste(wave_eqs, collapse = '\n')
  return(all_wave_model)
}

#' create_lav_ints
#'
#' @param manifest_vec
#' @param group_equal
#' @param long_free
#' @param num_groups
#' @param fix_first_intercept
#'
#' @return
#' @export
create_lav_ints <- function(manifest_vec, group_equal = F, long_free = F, num_groups = 4,
                            fix_first_intercept = F){
  if(long_free && group_equal){
    stop("Can't set groups equal and free longitudinal constraints.")
  }
  int_grp_eqs <- lapply(manifest_vec, function(man_var){
    int_label <- paste0('int_', which(manifest_vec == man_var))
    grp_labels <- rep(int_label, num_groups)
    if(!group_equal){
      grp_labels <- paste0(grp_labels, 'g', 1:num_groups)
    }
    grp_labels_collapsed <- paste(grp_labels, collapse = ', ')
    eq_rhs <- paste0('c(', grp_labels_collapsed, ')*1')
    eq_lhs <- paste0(letters[1:4], man_var)
    if(long_free && which(manifest_vec == man_var) > 1){
      eq <- paste0(eq_lhs, ' ~ ', '1')
    } else if(fix_first_intercept && which(manifest_vec == man_var) == 1){
      eq <- paste0(eq_lhs, ' ~ ', '0*1')
    } else {
      eq <- paste0(eq_lhs, ' ~ ', eq_rhs)
    }
    paste(eq, collapse = '\n')
  })
  int_all_eqs <- paste(int_grp_eqs, collapse = '\n\n')
  return(int_all_eqs)
}

#' create_lav_resid_var
#'
#' @param manifest_vec
#' @param group_equal
#' @param long_free
#' @param num_groups
#'
#' @return
#' @export
create_lav_resid_var <- function(manifest_vec, group_equal = F, long_free = F, num_groups = 4){
  if(long_free && group_equal){
    stop("Can't set groups equal and free longitudinal constraints.")
  }
  var_grp_eqs <- lapply(manifest_vec, function(man_var){
    var_label <- paste0('v_', which(manifest_vec == man_var))
    grp_labels <- rep(var_label, num_groups)
    if(!group_equal){
      grp_labels <- paste0(grp_labels, 'g', 1:num_groups)
    }
    grp_labels_collapsed <- paste(grp_labels, collapse = ', ')
    eq_rhs <- paste0('c(', grp_labels_collapsed, ')')
    eq_lhs <- paste0(letters[1:4], man_var)
    if(long_free){
      eq <- paste0(eq_lhs, ' ~~ ', eq_lhs)
    } else {
      eq <- paste0(eq_lhs, ' ~~ ', eq_rhs, '*', eq_lhs)
    }
    paste(eq, collapse = '\n')
  })
  var_all_eqs <- paste(var_grp_eqs, collapse = '\n\n')
  return(var_all_eqs)
}

#' create_lav_resid_covar
#'
#' @param manifest_vec
#' @param group_equal
#' @param long_free
#' @param num_groups
#'
#' @return
#' @export
create_lav_resid_covar <- function(manifest_vec, group_equal = F, long_free = F, num_groups = 4){
  if(long_free && group_equal){
    stop("Can't set groups equal and free longitudinal constraints.")
  }
  covar_grp_eqs <- lapply(manifest_vec, function(man_var){
    lag_distances <- unlist(lapply(3:1, function(x) seq(from = 1, to = x, by = 1)))
    var_label <- paste0('cv_', which(manifest_vec == man_var), lag_distances)
    grp_labels_collapsed <- lapply(var_label, function(alab){
      thelabs <- rep(alab, num_groups)
      if(!group_equal){
        thelabs <- paste0(thelabs, 'g', 1:num_groups)
      }
      thelabs_collapsed <- paste(thelabs, collapse = ', ')
      return(thelabs_collapsed)
    })
    if(long_free){
      connectors <- ' ~~ '
    } else {
      connectors <- paste0(' ~~ c(', grp_labels_collapsed, ')*')
    }
    wave_vars <- paste0(letters[1:4], man_var)
    covar_label_matrix <- combn(wave_vars, 2)
    eqs <- paste0(covar_label_matrix[1,], connectors, covar_label_matrix[2,])
    eqs_collapsed <- paste(eqs, collapse = '\n')
    return(eqs_collapsed)
  })
  covar_grp_eqs_all <- paste(covar_grp_eqs, collapse = '\n\n')
  return(covar_grp_eqs_all)
}

#' create_lav_invariance_model
#'
#' @param factor_name
#' @param manifest_vec
#' @param type
#' @param num_groups
#' @param no_widman
#' @param fix_first_intercept
#'
#' @return
#' @export
create_lav_invariance_model <- function(factor_name, manifest_vec, type = c('unconstrained', 'long_metric', 'long_strong', 'long_strict', 'baseline','metric','strong', 'strict'), num_groups = 4, no_widman = F, fix_first_intercept = F){
  #unconstrained: estimate first loading, constrain across waves, but freely
  #estimate remaining loadings. Estimate first intercept but constrain to be
  #invariant across time. Freely estimate remaining parameters. Allow covariance
  #among residuals of same item across waves.
  #
  #long_metric: constrain loadings to be the same across time (but vary across
  #group).
  #
  #long_strong: constain intercepts to be the same across time (but vary across
  #group).
  #
  #long_scrict: constain residuals to be the same across time within item (but
  #vary across group).
  #
  #baseline: strict longitudinal invariance with added constraint that
  #covariance of residuals across the same lag distance should be the same.
  #
  #metric: factor loadings equal across groups.
  #
  #strong: intercepts equal across groups.
  #
  #strict: residual variances are equal across groups.

  amod <- switch(type,
                 unconstrained = paste(
                   c(create_lav_factor_loadings(
                     factor_name, manifest_vec,
                     group_equal = F, long_free = T, num_groups = num_groups),
                     create_lav_ints(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups,
                       fix_first_intercept = fix_first_intercept),
                     create_lav_resid_var(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups),
                     create_lav_resid_covar(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups)),
                   collapse = '\n\n#---\n\n'),
                 long_metric = paste(
                   c(create_lav_factor_loadings(
                     factor_name, manifest_vec,
                     group_equal = F, long_free = F, num_groups = num_groups),
                     create_lav_ints(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups,
                       fix_first_intercept = fix_first_intercept),
                     create_lav_resid_var(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups),
                     create_lav_resid_covar(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups)),
                   collapse = '\n\n#---\n\n'),
                 long_strong = paste(
                   c(create_lav_factor_loadings(
                     factor_name, manifest_vec,
                     group_equal = F, long_free = F, num_groups = num_groups),
                     create_lav_ints(
                       manifest_vec,
                       group_equal = F, long_free = F, num_groups = num_groups,
                       fix_first_intercept = fix_first_intercept),
                     create_lav_resid_var(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups),
                     create_lav_resid_covar(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups)),
                   collapse = '\n\n#---\n\n'),
                 long_strict = paste(
                   c(create_lav_factor_loadings(
                     factor_name, manifest_vec,
                     group_equal = F, long_free = F, num_groups = num_groups),
                     create_lav_ints(
                       manifest_vec,
                       group_equal = F, long_free = F, num_groups = num_groups,
                       fix_first_intercept = fix_first_intercept),
                     create_lav_resid_var(
                       manifest_vec,
                       group_equal = F, long_free = F, num_groups = num_groups),
                     create_lav_resid_covar(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups)),
                   collapse = '\n\n#---\n\n'),
                 baseline = paste(
                   c(create_lav_factor_loadings(
                     factor_name, manifest_vec,
                     group_equal = F, long_free = F, num_groups = num_groups),
                     create_lav_ints(
                       manifest_vec,
                       group_equal = F, long_free = F, num_groups = num_groups,
                       fix_first_intercept = fix_first_intercept),
                     create_lav_resid_var(
                       manifest_vec,
                       group_equal = F, long_free = F, num_groups = num_groups),
                     create_lav_resid_covar(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups)),
                   collapse = '\n\n#---\n\n'),
                 metric = paste(
                   c(create_lav_factor_loadings(
                     factor_name, manifest_vec, group_equal = T, num_groups = num_groups),
                     create_lav_ints(
                       manifest_vec, group_equal = F, num_groups = num_groups,
                       fix_first_intercept = fix_first_intercept),
                     create_lav_resid_var(
                       manifest_vec, group_equal = F, num_groups = num_groups),
                     create_lav_resid_covar(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups)),
                   collapse = '\n\n#---\n\n'),
                 strong = paste(
                   c(create_lav_factor_loadings(
                     factor_name, manifest_vec, group_equal = T, num_groups = num_groups),
                     create_lav_ints(
                       manifest_vec, group_equal = T, num_groups = num_groups,
                       fix_first_intercept = fix_first_intercept),
                     create_lav_resid_var(
                       manifest_vec, group_equal = F, num_groups = num_groups),
                     create_lav_resid_covar(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups)),
                   collapse = '\n\n#---\n\n'),
                 strict = paste(
                   c(create_lav_factor_loadings(
                     factor_name, manifest_vec, group_equal = T, num_groups = num_groups),
                     create_lav_ints(
                       manifest_vec, group_equal = T, num_groups = num_groups,
                       fix_first_intercept = fix_first_intercept),
                     create_lav_resid_var(
                       manifest_vec, group_equal = T, num_groups = num_groups),
                     create_lav_resid_covar(
                       manifest_vec,
                       group_equal = F, long_free = T, num_groups = num_groups)),
                   collapse = '\n\n#---\n\n'))
  if(no_widman){
    return(amod)
  } else {
    amod <- paste(create_lav_id_strat_widaman(factor_name = factor_name), amod, sep = '\n')
    return(amod)
  }
}

#' get_factor_items
#'
#' @param scoring_df
#' @param factor_name
#'
#' @return
#' @export
get_factor_items <- function(scoring_df, factor_name){
  return(scoring_df$item[scoring_df$factor_name == factor_name])
}

#' get_item_colnames
#'
#' @param manifest_vec
#'
#' @return
#' @export
get_item_colnames <- function(manifest_vec){
  return(unlist(lapply(letters[1:4], paste0, manifest_vec)))
}

#' sem_invar
#'
#' @param model
#' @param data
#' @param ...
#'
#' @return
#' @export
#' @import lavaan
sem_invar <- function(model, data, ...){
  fit <- lavaan(model = model, data = data,
                int.lv.free = T,
                auto.fix.first = F,
                int.ov.free = TRUE,
                auto.var = TRUE,
                auto.cov.lv.x = TRUE,
                auto.cov.y = TRUE,
                ...)
  return(fit)
}

#' growth_cofs
#'
#' @param model
#' @param data
#' @param ...
#'
#' @return
#' @export
growth_cofs <- function(model, data, ...){
  fit <- lavaan::lavaan(model = model, data = data,
                        meanstructure = TRUE,
                        int.ov.free = TRUE,
                        int.lv.free = TRUE,
                        auto.fix.first = TRUE,
                        auto.fix.single = TRUE,
                        auto.var = TRUE,
                        auto.cov.lv.x = TRUE,
                        auto.th = TRUE,
                        auto.delta = TRUE,
                        auto.cov.y = TRUE, ...)
  return(fit)
}

#' run_invariance_model
#'
#' @param factor_name
#' @param manifest_vec
#' @param item_data
#' @param type
#' @param group
#' @param num_groups
#' @param ...
#'
#' @return
#' @export
run_invariance_model <- function(factor_name, manifest_vec, item_data, type = 'baseline', group = NULL, num_groups = 4, ...){
  an_invar_model <- create_lav_invariance_model(
    factor_name,
    manifest_vec,
    type = type,
    num_groups = num_groups)

  if(is.null(group)){
    a_fit <- try(sem_invar(model = an_invar_model,
                           data = item_data,
                           ...))
  } else {
    a_fit <- try(sem_invar(model = an_invar_model,
                           data = item_data,
                           group = group,
                           ...))
  }
  return(a_fit)
}

#' run_invariance_tests
#'
#' @param factor_name
#' @param manifest_vec
#' @param item_data
#' @param group
#' @param ...
#'
#' @return
#' @export
#' @import future
#' @import lavaan
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @import listenv
run_invariance_tests <- function(factor_name, manifest_vec, item_data, group, ...){
  all_invar_types <- c(
    'unconstrained',
    'long_metric',
    'long_strong',
    'long_strict',
    # 'baseline',
    'metric',
    'strong',
    'strict')

  invar_models_future_listenv <- listenv::listenv()
  for(it in seq_along(all_invar_types)){
    invar_models_future_listenv[[it]] <- future::future({
      message(paste0("Running: ", factor_name, ", ", all_invar_types[[it]]))
      try(run_invariance_model(factor_name,
                               manifest_vec,
                               item_data,
                               type = all_invar_types[[it]],
                               group = group,
                               ...))
    })
  }
  return(invar_models_future_listenv)
}

#' process_invariance_tests
#'
#' @param invar_models_list
#' @param factor_name
#' @param fit_measures
#'
#' @return
#' @export
process_invariance_tests <- function(invar_models_list, factor_name, fit_measures){
  all_invar_types <- c(
    'unconstrained',
    'long_metric',
    'long_strong',
    'long_strict',
    # 'baseline',
    'metric',
    'strong',
    'strict')
  message(paste0('Processing: ', factor_name))
  fit_statistics <- future.apply::future_lapply(
    invar_models_list,
    function(x){
      tryfitmeasure <- try(fitmeasures(x, fit.measures = fit_measures))

      if(any(class(tryfitmeasure) == 'try-error')){
        rez <- rep(NA, length(fit_measures))
        names(rez) <- fit_measures
        return(rez)
      } else {
        return(tryfitmeasure)
      }
    }
  )
  fit_statistics_len <- lapply(fit_statistics, length)
  fit_stat_df <- data_frame(
    factor_name = factor_name,
    invar_type = rep(all_invar_types, fit_statistics_len),
    fit_stat = names(unlist(fit_statistics)),
    value = unlist(fit_statistics))

  converged_index <- unlist(lapply(invar_models_list, lavaan::inspect, what = 'converged'))
  lrt_test <- do.call(
    function(...){
      try(lavTestLRT(..., model.names = all_invar_types[converged_index]))
    },
    invar_models_list[converged_index])

  if(any(class(lrt_test) == 'try-error') || !(any(class(lrt_test) == 'anova'))){
    lrt_test_df <- data.frame()
  } else {
    lrt_test_df <- tibble::as_data_frame(lrt_test)
    lrt_test_df$invar_type <- rownames(lrt_test)
    lrt_test_df <- tidyr::gather(lrt_test_df, fit_stat, value,
                                 -invar_type)
    lrt_test_df$factor_name <- factor_name
  }

  summary_df <- dplyr::bind_rows(fit_stat_df, lrt_test_df)

  names(invar_models_list) <- all_invar_types
  return(summary_df)
}

#' test_factors_for_invariance
#'
#' @param factor_names
#' @param indicators_df
#' @param item_data
#' @param fit_measures
#' @param group
#' @param save_fits
#' @param ...
#'
#' @return
#' @export
test_factors_for_invariance <- function(factor_names, indicators_df, item_data, fit_measures = c('mfi', 'cfi', 'rmsea', 'rmsea.ci.lower', 'rmsea.ci.upper'), group, save_fits = F, ...){
  require(future)
  invar_mods_future_listenv <- listenv::listenv()
  for(fn in seq_along(factor_names)){
    factor_name <- factor_names[[fn]]
    factor_items <- get_factor_items(indicators_df, factor_name)
    if(!check_manifests_exist(factor_items, item_data)) {
      stop(factor_name, ' does not have item columns in item_data')
    }
    this_factor_data <- item_data[, c(get_item_colnames(factor_items), group)]
    invar_mods_future_listenv[[fn]] <- run_invariance_tests(
      factor_name = factor_name,
      manifest_vec = factor_items,
      item_data = this_factor_data,
      group = group,
      ...)
  }
  procd_invar_stats_listenv <- listenv::listenv()
  for(fn in seq_along(factor_names)){
    factor_name <- factor_names[[fn]]
    a_factor_invariance_suite <- lapply(as.list(invar_mods_future_listenv[[fn]]), future::value)
    if(save_fits){
      message(paste0('Saving: ', factor_name))
      saveRDS(a_factor_invariance_suite,
              file.path(data_dir, paste0(factor_name, '_invariance_mods.RDS')))
    }
    procd_invar_stats_listenv[[fn]] <- process_invariance_tests(
      invar_models_list = a_factor_invariance_suite,
      factor_name = factor_name,
      fit_measures = fit_measures)
  }
  return(as.list(procd_invar_stats_listenv))
}
