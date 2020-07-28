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

#' fit.jagsMANEC
#'
#' Fits a variety of NEC models using jags and provides a model averaged predictions based on DIC model weights
#'
#' @inheritParams fit.jagsNEC
#'
#' @param model.set A vector of the names of model types to be fit. Currently defaults to
#' all available model types. If "NEC" is supplied, only the NEC models will be fit. If "ECx" is supplied,
#'  only continuous curve models will be fit. see ?fit.jagsNEC to see the available models to fit.
#'
#' @export
#' @return All successully fitted jags model fits, mod.stats a data.frame of model fit statistics, NEC a model
#' averaged posterior of the estimated NEC, and pred.vals a list of model averaged predictions.

fit.jagsMANEC <- function(data,
                          x.var,
                          y.var,
                          trials.var = NA,
                          x.type = NA,
                          y.type = NA,
                          burnin = 5000,
                          n.iter = burnin + 500,
                          n.iter.update = 10000,
                          n.tries = 3,
                          params = c("top", "beta", "NEC", "SS", "SSsim"),
                          over.disp = FALSE,
                          sig.val = 0.01,
                          model.set = "all",
                          ...) {
  check_inputs(data = data, x.var = x.var, y.var = y.var, trials.var = trials.var)

  if (model.set[1] == "NEC") {
    model.set <- c("NEC3param", "NEC4param", "NECHormesis", "NECsigmoidal")
  }
  if (model.set[1] == "ECx") {
    model.set <- c("ECx4param", "ECxWeibull1", "ECxWeibull2", "ECxLinear")
  }
  if (model.set[1] == "all") {
    model.set <- c(
      "NEC3param", "NEC4param", "NECHormesis", "NECsigmoidal",
      "ECxLinear", "ECxExp", "ECxsigmoidal",
      "ECx4param", "ECxWeibull1", "ECxWeibull2"
    )
  }
  if (model.set[1] == "bot_free") {
    model.set <- c(
      "NEC3param", "NECHormesis", "NECsigmoidal",
      "ECxLinear", "ECxExp", "ECxsigmoidal"
    )
  }

  # Fit each of the models
  mod.fits <- vector(mode = "list", length = length(model.set))
  names(mod.fits) <- model.set

  for (m in seq_along(model.set)) {
    model <- model.set[m]
    fit.m <- try(
      fit.jagsNEC(
        data = data,
        x.var = x.var,
        y.var = y.var,
        trials.var = trials.var,
        x.type = x.type,
        y.type = y.type,
        burnin = burnin,
        n.iter = n.iter,
        n.iter.update = n.iter.update,
        n.tries = n.tries,
        params = params,
        over.disp = over.disp,
        model = model
      ),
      silent = TRUE
    )
    if (!inherits(fit.m, "try-error")) {
      mod.fits[[m]] <- fit.m
    } else {
      mod.fits[[m]] <- NA
    }
  }



  # collate all the elements
  export.list <- c(
    extract_modStats(mod.fits),
    list(data = data, x.var = x.var, y.var = y.var, trials.var = trials.var, over.disp = over.disp)
  )
  # assign a class to the output
  class(export.list) <- "jagsMANECfit"

  # return the collated output
  return(export.list)
}
