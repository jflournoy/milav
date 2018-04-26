## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = T, message = F, error = T, echo = F,
  collapse = T,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(milav)
library(future)
library(semTools)
if(grepl('(^n\\d|talapas-ln1)', system('hostname', intern = T))){
  nworkers <- 28
  data_dir <- '/gpfs/projects/dsnlab/flournoy/data/mi_sim_results'
  niter <- 1000
  plan(tweak(multiprocess, gc = T, workers = nworkers))
} else {
  nworkers <- future::availableCores()
  data_dir <- '/data/jflournoy/simrez/'
  niter <- 2
  plan(tweak(multiprocess, gc = T, workers = nworkers))
}
message('Workers: ', nworkers)

if(!dir.exists(data_dir)){
  dir.create(data_dir)
}

mitest <- function(generative_mod, unconstrained_mod, constrained_mod, fit_measures, more_fit_measures = NULL, sample_nobs = c(2.5e2,2.5e2,2.5e2,2.5e2), return_all = F){
  simdata <- lavaan::simulateData(generative_mod, 
                                  sample.nobs = sample_nobs,
                                  int.lv.free = T, 
                                  auto.fix.first = F, int.ov.free = TRUE, auto.var = TRUE, 
                                  auto.cov.lv.x = TRUE, auto.cov.y = TRUE,
                                  standardized = TRUE, empirical = F)
  unconstrained_fit <- milav::sem_invar(model = unconstrained_mod, data = simdata, group = 'group')
  constrained_fit <- milav::sem_invar(model = constrained_mod, data = simdata, group = 'group')
  fit_diff <- diff(
    do.call(rbind, 
            lapply(list(unconstrained_fit, constrained_fit),
                   lavaan::fitmeasures, fit.measures = fit_measures)))
  if(!is.null(more_fit_measures)){
    more_fit_diff <- diff(
    do.call(rbind, 
            lapply(list(unconstrained_fit, constrained_fit),
                   semTools::moreFitIndices, fit.measures = more_fit_measures)))
    fit_diff <- cbind(fit_diff, more_fit_diff)
  }
  if(return_all){
    return(list(fit_diff = fit_diff, 
                data = simdata,
                generative_mod = generative_mod,
                unconstrained_fit = unconstrained_fit, 
                constrained_fit = constrained_fit))
  }
  else{
    return(fit_diff)
  }
}

create_generative_mod <- function(nitems, loadings_from, loadings_to){
  loadings <- as.character(round(runif(nitems, loadings_from, loadings_to),2))
  names(loadings) <- paste0('c\\((L',1:nitems,',* *)+\\)')
  template <- milav::create_lav_invariance_model(factor_name = 'F',
                                                 manifest_vec = paste0('x', 1:nitems),
                                                 'strict')
  generative_mod <- stringr::str_replace_all(template, loadings)
  return(generative_mod)
}

run_mi_sims <- function(niter, nitems, loadings_from, loadings_to, unconstrained_mod, constrained_mod, fit_measures, more_fit_measures = NULL, sample_nobs = c(2.5e2,2.5e2,2.5e2,2.5e2), return_all = F) {
  require(future)
  require(listenv)
  fit_diffs <- listenv::listenv()
  for(i in 1:niter){
    fit_diffs[[i]] <- future({
      message(paste(c('Iter: ', ', L_f: ', ', L_t: ', ', N Items: '), 
                    c(i, loadings_from, loadings_to, nitems)))
      generative_mod <- create_generative_mod(nitems = nitems,
                                              loadings_from = loadings_from,
                                              loadings_to = loadings_to)
      
      fit_diff <- mitest(generative_mod = generative_mod,
                         unconstrained_mod = unconstrained_mod,
                         constrained_mod = constrained_mod,
                         fit_measures = fit_measures,
                         more_fit_measures = more_fit_measures,
                         sample_nobs = sample_nobs, 
                         return_all = return_all)
      fit_diff
    })   
  }
  return(fit_diffs)
}

fit_measures <- c('mfi', 'cfi', 'rmsea')
more_fit_measures <- c('gammaHat', 'adjGammaHat')
savefile <- file.path(data_dir, 'mi_simulation_results.RDS')

properties_df <- expand.grid(loadings_from = seq(.5, .9, .1), nitems = c(4,6,8))
properties_df$loadings_to <- properties_df$loadings_from + .09

if(!file.exists(savefile)){
  l_fut <- listenv::listenv()
  for(i in 1:nrow(properties_df)){
    long_strict_mod <- milav::create_lav_invariance_model(
      factor_name = 'F',
      manifest_vec = paste0('x', 1:properties_df$nitems[[i]]),
      'long_strict')
    group_strict_mod <- milav::create_lav_invariance_model(
      factor_name = 'F',
      manifest_vec = paste0('x', 1:properties_df$nitems[[i]]),
      'strict')
    l_fut[[i]] <- run_mi_sims(
      niter = niter, nitems = properties_df$nitems[[i]], 
      loadings_from = properties_df$loadings_from[[i]], 
      loadings_to = properties_df$loadings_to[[i]], 
      unconstrained_mod = long_strict_mod, constrained_mod = group_strict_mod,
      fit_measures = fit_measures,
      more_fit_measures = more_fit_measures)
  }
  
  l_rez <- lapply(l_fut, lapply, future::value)
  names(l_rez) <- paste0('Loadings: ', 
                         properties_df$loadings_from, ' - ', 
                         properties_df$loadings_to,
                         '; N items: ', 
                         properties_df$nitems)
  saveRDS(l_rez, file = file.path(data_dir, 'mi_simulation_results.RDS'))
  
} else {
  l_rez <- readRDS(savefile)
}

## ----eval = T, fig.width=6, fig.height=5---------------------------------
long_strict_mod <- milav::create_lav_invariance_model(
  factor_name = 'F',
  manifest_vec = paste0('x', 1:properties_df$nitems[[1]]),
  'long_strict')
group_strict_mod <- milav::create_lav_invariance_model(
  factor_name = 'F',
  manifest_vec = paste0('x', 1:properties_df$nitems[[1]]),
  'strict')
asim_f <- run_mi_sims(
      niter = niter, nitems = properties_df$nitems[[1]], 
      loadings_from = properties_df$loadings_from[[1]], 
      loadings_to = properties_df$loadings_to[[1]], 
      unconstrained_mod = long_strict_mod, constrained_mod = group_strict_mod,
      fit_measures = fit_measures,
      more_fit_measures = more_fit_measures, return_all = T)
asim <- value(asim_f[[1]])
semPlot::semPaths(asim$constrained_fit, include = 1, ask = F)

## ----fig.width=6, fig.height=5-------------------------------------------
all_df <- tidyr::extract(
  dplyr::bind_rows(lapply(l_rez, function(l) as.data.frame(do.call(what = rbind, t(l)))),
                           .id = 'model'),
  model, 
  c('loading', 'nitems'), 
  'Loadings: (.*); N items: (\\d+)')

all_df_quant <- dplyr::summarise_at(
  dplyr::group_by(all_df, loading, nitems),
  dplyr::vars('gammaHat', 'mfi', 'cfi'),
  dplyr::funs(
    '.005' = quantile(., probs = c(.005)),
    '.01' = quantile(., probs = c(.01)),
    '.05' = quantile(., probs = c(.05)),
    '.50' = quantile(., probs = c(.50))))

ggplot2::ggplot(all_df_quant,
                ggplot2::aes(x = as.numeric(nitems), y = cfi_.01, linetype = loading))+
  ggplot2::geom_ribbon(ggplot2::aes(ymin = cfi_.005, ymax = cfi_.05),
                       alpha = .2)+
  ggplot2::geom_point()+
  ggplot2::geom_line()+
  ggplot2::geom_text(ggplot2::aes(x = 8.3, 
                                  y = cfi_.01, 
                                  label = loading), 
                     data = all_df_quant[all_df_quant$nitems==8,],
                     size = 3.5)+
  ggplot2::coord_cartesian(xlim = c(4, 8.5)) +
  ggplot2::theme(legend.position="none") + 
  ggplot2::labs(title = 'CFI distribution under null of invariance',
                subtitle = '1000 iterations per condition',
                x = 'Number of items',
                y = 'CFI',
                caption = 'Lines are at 1%-ile, labeled by range of loadings. 
                Shaded regions encompass the distribution range of [5%, 0.5%]')


ggplot2::theme_set(ggplot2::theme_minimal())
ggplot2::ggplot(all_df_quant,
                ggplot2::aes(x = as.numeric(nitems), y = mfi_.01, linetype = loading))+
  ggplot2::geom_ribbon(ggplot2::aes(ymin = mfi_.005, ymax = mfi_.05),
                       alpha = .2)+
  ggplot2::geom_point()+
  ggplot2::geom_line()+
  ggplot2::geom_text(ggplot2::aes(x = 8.3, 
                                  y = mfi_.01+c(.00032, -.00029, -.0000, 0, 0), 
                                  label = loading), 
                     data = all_df_quant[all_df_quant$nitems==8,],
                     size = 3.5)+
  ggplot2::coord_cartesian(ylim = c(0,-.0175), xlim = c(4, 8.5)) +
  ggplot2::scale_y_continuous(breaks = c(0, -.007, -.01, -.012, -.014)) + 
  ggplot2::theme(legend.position="none") + 
  ggplot2::labs(title = 'MFI distribution under null of invariance',
                subtitle = '1000 iterations per condition',
                x = 'Number of items',
                y = 'McDonald\'s noncentrality index',
                caption = 'Lines are at 1%-ile, labeled by range of loadings. 
                Shaded regions encompass the distribution range of [5%, 0.5%]')


ggplot2::ggplot(all_df_quant,
                ggplot2::aes(x = as.numeric(nitems), y = gammaHat_.01, linetype = loading))+
  ggplot2::geom_ribbon(ggplot2::aes(ymin = gammaHat_.005, ymax = gammaHat_.05),
                       alpha = .2)+
  ggplot2::geom_point()+
  ggplot2::geom_line()+
  ggplot2::geom_text(ggplot2::aes(x = 8.3, 
                                  y = gammaHat_.01+c(.00009, -.00000, -.0000, 0, 0), 
                                  label = loading), 
                     data = all_df_quant[all_df_quant$nitems==8,],
                     size = 3.5)+
  ggplot2::coord_cartesian(ylim = c(0,-.0045), xlim = c(4, 8.5)) +
  ggplot2::scale_y_continuous(breaks = c(0, -.0015, -.0025, -.003, -.0004)) + 
  ggplot2::theme(legend.position="none") + 
  ggplot2::labs(title = expression(paste('Distribution of ', hat(gamma), ' under null of invariance')),
                subtitle = '1000 iterations per condition',
                x = 'Number of items',
                y = expression(hat(gamma)),
                caption = 'Lines are at 1%-ile, labeled by range of loadings. 
                Shaded regions encompass the distribution range of [5%, 0.5%]')



## ----results='asis'------------------------------------------------------
all <- do.call(rbind, lapply(l_rez, do.call, what = rbind))

knitr::kable(apply(all, 2, quantile, probs = c(.005/(8), .005, .01, .05, .5, .95, .99)),
             caption = 'Quantiles over all iterations')

for(i in seq_along(l_rez)){
  cat('\n\n')
  cat(paste(knitr::kable(
    apply(do.call(rbind, l_rez[[i]]), 
          2, quantile, 
          probs = c(.01, .05, .5, .95, .99)),
    caption = names(l_rez)[[i]], digits = 4), collapse = '\n'))
}

## ----echo=F--------------------------------------------------------------
plot_mi_rez <- function(df, xvar, i, n){
  mfi.q <- quantile(df[,which(colnames(df) == xvar)], probs = c(.01, .05, .5))
  aplot <- ggplot2::ggplot(
    as.data.frame(df),
    ggplot2::aes_string(x = xvar))+
      # ggplot2::geom_histogram(binwidth = .0005) + 
      ggplot2::geom_density() + 
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = `1%`, linetype = '1%'), 
        data = as.data.frame(t(mfi.q))) +
      ggplot2::geom_text(
        ggplot2::aes(x = `1%`, y = 50, label = round(`1%`, 3)),
        nudge_x = -.0015, alpha = .8,
        data = as.data.frame(t(mfi.q))) +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = `5%`, linetype = '5%'), 
        data = as.data.frame(t(mfi.q))) +
      ggplot2::geom_text(
        ggplot2::aes(x = `5%`, y = 50, label = round(`5%`, 3)),
        nudge_x = -.0015, alpha = .8,
        data = as.data.frame(t(mfi.q))) +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = `50%`, linetype = '50%'), 
        data = as.data.frame(t(mfi.q))) +
      ggplot2::geom_text(
        ggplot2::aes(x = `50%`, y = 50, label = round(`50%`, 3)),
        nudge_x = -.0015, alpha = .8,
        data = as.data.frame(t(mfi.q))) +
      ggplot2::scale_linetype_manual(values = c('1%' = 2, '5%' = 3, '50%' = 4),
                                     name = 'Percentile') + 
      ggplot2::theme_minimal() +
      ggplot2::labs(title = n[[i]])
  return(aplot)
}

## ----fig.width=6, fig.height=5-------------------------------------------
plot_mi_rez(all, xvar = 'mfi', i = 1, c('All MFI'))
plot_mi_rez(all, xvar = 'gammaHat', i = 1, c('All gammaHat'))

## ----mfiplots, fig.width=6, fig.height=5---------------------------------
nope <- lapply(seq_along(l_rez), function(i, l, n){
  print(plot_mi_rez(do.call(rbind, l[[i]]), 'mfi', i, n))
}, l = l_rez, n = names(l_rez))

## ----cfiplots, fig.width=6, fig.height=5---------------------------------
nope <- lapply(seq_along(l_rez), function(i, l, n){
  print(plot_mi_rez(do.call(rbind, l[[i]]), 'cfi', i, n))
}, l = l_rez, n = names(l_rez))

## ----gammahatplots, fig.width=6, fig.height=5----------------------------
nope <- lapply(seq_along(l_rez), function(i, l, n){
  print(plot_mi_rez(do.call(rbind, l[[i]]), 'gammaHat', i, n))
}, l = l_rez, n = names(l_rez))

## ----adjgammahatplots, fig.width=6, fig.height=5-------------------------
nope <- lapply(seq_along(l_rez), function(i, l, n){
  print(plot_mi_rez(do.call(rbind, l[[i]]), 'adjGammaHat', i, n))
}, l = l_rez, n = names(l_rez))

