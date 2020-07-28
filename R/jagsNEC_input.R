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

#' jagsNEC_input
#'
#' Checks the model inputs for a jagsNEC model fit
#'
#' @inheritParams fit.jagsNEC
#'
#' @details
#'
#' This is a wrapper function to test input data criteria and write the jags model file for use in a jagsNEC model fit
#'
#' @export
#' @return Modified elements of the jagsNEC input data.

jagsNEC_input <- function(data,
                          x.var,
                          y.var,
                          trials.var,
                          x.type,
                          y.type,
                          params,
                          over.disp,
                          model) {
  if (is.na(y.type) == F) {
    if (over.disp == TRUE & y.type == "beta") {
      y.type <- NA
    }
  }

  check_inputs(data = data, x.var = x.var, y.var = y.var, trials.var = trials.var)
  # extract the data
  y.dat <- data[, y.var]
  x.dat <- data[, x.var]


  # check the data are lower at high x compared to low x (ie the response variable declines within increase in the x)
  if (mean(y.dat[which(x.dat < mean(x.dat))]) < mean(y.dat[which(x.dat > mean(x.dat))]) & model != "NECHormesis") {
    stop("The mean value of the response for the lower half of the
             concentration data are lower than that of the upper half of the concentration data.
             jagsNEC only fits concentration response data where the
             response declines with increasing values of concentration.")
  }

  # check variable type x.var
  if (is.na(x.type) == TRUE) { # if the x.var is not specified, then guess
    if (class(x.dat) == "integer") {
      stop("jagsNEC does not currently support integer concentration data. Please provide
           a numeric x.var")
    }
    if (class(x.dat) == "numeric" & max(x.dat) > 1 & min(x.dat) >= 0) {
      x.type <- "gamma"
    }
    if (class(x.dat) == "numeric" & max(x.dat) <= 1 & min(x.dat) >= 0) {
      x.type <- "beta"
    }
    if (class(x.dat) == "numeric" & min(x.dat) < 0) {
      x.type <- "gaussian"
    }
  }

  # check variable type y.var
  if (is.na(y.type) == T) { # if the y.var is not specified, then guess
    if (class(y.dat) == "numeric" & max(y.dat) > 1 & min(y.dat) >= 0) {
      y.type <- "gamma"
    }
    if (class(y.dat) == "numeric" & max(y.dat) <= 1 & min(y.dat) >= 0) {
      y.type <- "beta"
    }
    if (class(y.dat) == "numeric" & min(y.dat) < 0) {
      y.type <- "gaussian"
    }
    if (class(y.dat) == "integer" & min(y.dat) >= 0 & is.na(trials.var) == TRUE) {
      y.type <- "poisson"
    }
    if (is.na(trials.var) != TRUE & class(y.dat) != "integer") {
      stop("You have supplied a trials.var argument, suggesting you wish to model a binomial.
             Please ensure y.var is an integer representing the number of successes.")
    }

    if (class(y.dat) == "integer" & min(y.dat) >= 0 & is.na(trials.var) != TRUE) {
      y.type <- "binomial"
    }
  }

  # check there is a valid model type
  if (is.na(match(model, c(
    "NEC3param", "NECsigmoidal", "NEC4param", "NECHormesis",
    "ECx4param", "ECxWeibull1", "ECxWeibull2", "ECxLinear", "ECxExp", "ECxsigmoidal"
  )))) {
    stop("The model type you have specified does not extist.")
  }

  if (y.type == "poisson" & over.disp == TRUE) {
    y.type <- "negbin"
  }
  if (y.type == "binomial" & over.disp == TRUE) {
    y.type <- "beta"
    data[, y.var] <- data[, y.var] / data[, trials.var]
  }

  if (y.type == "gamma") {
    params <- c(params, "shape")
  }
  if (y.type == "gaussian") {
    params <- c(params, "alpha", "sigma")
  }
  if (y.type == "negbin") {
    params <- c(params, "size")
  }

  # error catching for 0 for gamma by adding very small value (no tweedie available in jags)
  if (min(data[, x.var]) == 0 & x.type == "gamma") {
    tt <- data[, x.var]
    min.val <- min(tt[which(tt > 0)])
    data[which(tt == 0), x.var] <- tt[which(tt == 0)] + (min.val / 100)
  }

  if (min(data[, y.var]) == 0 & y.type == "gamma") {
    tt <- data[, y.var]
    min.val <- min(tt[which(tt > 0)])
    data[which(tt == 0), y.var] <- tt[which(tt == 0)] + (min.val / 100)
  }
  # error catching for 0 for beta by adding very small value (beta does not take zero)
  if (min(data[, x.var]) == 0 & x.type == "beta") {
    tt <- data[, x.var]
    min.val <- min(tt[which(tt > 0)])
    data[which(tt == 0), x.var] <- tt[which(tt == 0)] + (min.val / 100)
  }

  if (min(data[, y.var]) == 0 & y.type == "beta") {
    tt <- data[, y.var]
    min.val <- min(tt[which(tt > 0)])
    data[which(tt == 0), y.var] <- tt[which(tt == 0)] + (min.val / 100)
  }

  # error catching for 1 for beta by subtracting very small value (beta does not take 1)
  if (max(data[, x.var]) == 1 & x.type == "beta") {
    tt <- data[, x.var]
    data[which(tt == 1), x.var] <- tt[which(tt == 1)] - 0.001
  }

  if (max(data[, y.var]) == 1 & y.type == "beta") {
    tt <- data[, y.var]
    data[which(tt == 1), y.var] <- tt[which(tt == 1)] - 0.001
  }

  # create jags model data list
  mod.dat <<- list(
    x = data[, x.var], # concentration
    y = data[, y.var], # response (successes)

    N = nrow(data)
  ) # Sample size


  response <- data[, y.var]

  if (y.type == "binomial") {
    mod.dat$trials <- data[, trials.var] # number of "trials"
    response <- data[, y.var] / data[, trials.var]
    if (max(response) > 1) {
      stop(paste("Your successes as indicated in ", y.var, " exceed the number of trials contained in ", trials.var, ".", sep = ""))
    }
  }

  mod.file <- Write.jagsModFile(x.type, y.type, mod.dat, params, model)
  init.fun <- mod.file$init.fun
  params <- mod.file$params

  return(list(
    params = params,
    response = response,
    mod.dat = mod.dat,
    data = data,
    y.type = y.type,
    x.type = x.type,
    x.dat = x.dat,
    y.dat = y.dat,
    init.fun = init.fun
  ))
}
