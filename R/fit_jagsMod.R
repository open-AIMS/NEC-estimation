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

#' fit.jagsMod
#'
#' Fits a Bayesian model using jags and R2jags
#'
#' @param mod.dat the model data to use for the model fit
#'
#' @param init.fun the init function to use for the model fit
#'
#' @inheritParams fit.jagsNEC
#'
#' @export
#' @return The fitted jags model output

fit.jagsMod <- function(mod.dat, init.fun, params, burnin, n.iter, init.value.warning, n.tries) {
  all.Js <- list()

  warn <- getOption("warn")
  options(warn = -1)
  J1 <- try(R2jags::jags(
    data = mod.dat,
    inits = init.fun,
    parameters = params,
    model = "NECmod.txt",
    n.thin = 10,
    n.chains = 5,
    n.burnin = burnin,
    n.iter = n.iter,
    progress.bar = "none"
  ), silent = T)
  options(warn = warn)

  if (class(J1) != "try-error") { # make sure the fitted model had good mixing
    all.Js <- c(all.Js, list(J1))
    if (max(unlist(check.mixing(J1)$cv.test)) == 1) {
      class(J1) <- "try-error"
    }
  }
  if (init.value.warning == TRUE) {
    warning("The generated init.fun failed to yield a valid model. Model based on default jags initial values")
  }
  # if the attempt fails try 10 more times
  warn <- getOption("warn")
  options(warn = -1)
  w <- 1
  while (class(J1) == "try-error" & w <= n.tries) {
    w <- w + 1
    J1 <- try(R2jags::jags(
      data = mod.dat,
      inits = init.fun,
      parameters = params,
      model = "NECmod.txt",
      n.thin = 10,
      n.chains = 5,
      n.burnin = burnin,
      n.iter = n.iter,
      progress.bar = "none"
    ), silent = T)
    if (class(J1) != "try-error") { # make sure the fitted model had good mixing
      all.Js <- c(all.Js, list(J1))
      if (max(unlist(check.mixing(J1)$cv.test)) == 1) {
        class(J1) <- "try-error"
      }
    }
  }

  # if the attempt fails try without initial values
  w <- 1
  while (class(J1) == "try-error" & w <= n.tries) {
    w <- w + 1
    J1 <- try(R2jags::jags(
      data = mod.dat,
      parameters = params,
      model = "NECmod.txt",
      n.thin = 10,
      n.chains = 5,
      n.burnin = burnin,
      n.iter = n.iter,
      progress.bar = "none"
    ), silent = T)
    if (class(J1) != "try-error") { # make sure the fitted model had good mixing
      all.Js <- c(all.Js, list(J1))
      if (max(unlist(check.mixing(J1)$cv.test)) == 1) {
        class(J1) <- "try-error"
      }
    }
  }
  # options(warn=warn)
  return(list(J1 = J1, all.Js = all.Js))
}
