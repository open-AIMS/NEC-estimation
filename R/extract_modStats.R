#    Copyright 2020 Australian Institute of Marine Science
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

#' extract.modStats
#'
#' Extracts a range of statistics from a list of jagsNECfit model fits.
#'
#' @param  mod.fits a jagsMANECfit mod.fits output list, as returned by fit.jagsMANEC
#'
#' @export
#' @return A list of model statistical output derived from the input model list

extract_modStats <- function(mod.fits) {
  model.set <- names(mod.fits)
  # extract model parameters that do not vary across models
  success.models <- model.set[sapply(mod.fits, FUN = class) == "jagsNECfit"]
  if (length(success.models) == 0) {
    stop("None of the models fit successfully, 
     try using fit.jagsNEC instead using the default settings as a starting point for trouble shooting.")
  } else {
    warning(paste("successfully fitted the models: ", paste(success.models, collapse = " ")))
  }

  mod.fits <- mod.fits[success.models]
  # extract the model statistics for each fit
  mod.stats <- data.frame(DIC = sapply(mod.fits, FUN = function(x) {
    x$DIC
  }))
  mod.stats$DIC.delta <- mod.stats$DIC - min(mod.stats$DIC, na.rm = TRUE)

  # sapply(mod.fits, FUN=function(x){loo::waic(x$BUGSoutput$sims.list$loglik)})


  # mod.stats$WAIC <-
  # mod.stats$WAIC.delta <- mod.stats$WAIC-min(mod.stats$WAIC, na.rm = TRUE)
  mod.stats$wi <- wi(mod.stats$DIC)
  mod.stats$pD <- unlist(lapply(mod.fits, FUN = function(x) {
    x$pD
  }))
  mod.stats$over.disp <- unlist(lapply(mod.fits, FUN = function(x) {
    x$over.disp
  }))

  mcmc.stats <- unique(do.call("rbind", lapply(mod.fits, FUN = function(x) {
    x[c("n.chains", "n.iter", "n.burnin", "n.thin", "n.keep", "n.sims", "y.type", "x.type")]
  })))
  mcmc.list <- do.call("list", mcmc.stats)
  names(mcmc.list) <- colnames(mcmc.stats)
  sample.size <- mcmc.list$n.sims

  # model averaged NEC posterior
  NEC.posterior <- unlist(lapply(success.models, FUN = function(x) {
    base::sample(mod.fits[[x]]$sims.list$NEC, round(mcmc.list$n.sims * mod.stats[x, "wi"]))
  }))

  # model averaged predicted y
  predicted.y <- rowSums(do.call("cbind", lapply(success.models, FUN = function(x) {
    mod.fits[[x]]$predicted.y * mod.stats[x, "wi"]
  })))

  # model averaged pred.vals
  x <- mod.fits[[success.models[1]]]$pred.vals$x

  y.m <- rowSums(do.call("cbind", lapply(success.models, FUN = function(x) {
    mod.fits[[x]]$pred.vals$y.m * mod.stats[x, "wi"]
  })))

  # model weighted posterior
  posterior.predicted <- do.call("cbind", lapply(success.models, FUN = function(x) {
    mod.fits[[x]]$pred.vals$posterior[
      ,
      base::sample(1:sample.size, round(sample.size * mod.stats[x, "wi"]))
    ]
  }))

  y <- apply(posterior.predicted, MARGIN = 1, FUN = median)
  up <- apply(posterior.predicted, MARGIN = 1, FUN = quantile, probs = 0.975)
  lw <- apply(posterior.predicted, MARGIN = 1, FUN = quantile, probs = 0.025)

  # collate all the elements
  export.list <-
    c(
      mcmc.list,
      list(
        mod.fits = mod.fits,
        success.models = success.models,
        mod.dat = mod.fits[[1]]$mod.dat,
        mod.stats = mod.stats,
        mcmc.stats = mcmc.stats,
        sims.list = list(NEC = NEC.posterior),
        predicted.y = predicted.y,
        residuals = mod.fits[[1]]$mod.dat$y - predicted.y,
        pred.vals = list(x = x, y = y, up = up, lw = lw, posterior = posterior.predicted, y.m = y.m),
        NEC = quantile(NEC.posterior, probs = c(0.025, 0.5, 0.975))
      )
    )
  return(export.list)
}
