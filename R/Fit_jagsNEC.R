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

#' fit.jagsNEC
#'
#' Fits an NEC model as per Fox 2010, using jags and R2jags
#'
#' @param  data a data.frame containing the data to use for the model
#'
#' @param x.var the column heading indicating the concentration (x) variable
#'
#' @param y.var the column heading indicating the response (y) variable
#'
#' @param trials.var the column heading indicating the column for the number of "trials" for binomial response data.
#' If not supplied, the model may run but will not be the model you intended!
#'
#' @param x.type the statistical distribution to use for the x (concentration) data. This will be guess based on the
#' characteristic of the input data if not supplied.
#'
#' @param y.type the statistical distribution to use for the y (response) data. This may currently be one of  'binomial',
#' 'poisson',' 'gaussian', or 'gamma'. Others can be added as required, please contact the package maintainer.
#' If not supplied, the appropriate distribution will be guessed based on the distribution of the input data.
#'
#' @param  params A vector of names indicating the parameters that to trace during the jags fit. For the NEC jags model
#' this is typically 'NEC','top' and 'beta'. If left out, fit.jagsNEC will supply this based on the selected y.type and
#' x.type.
#'
#' @param burnin the number of iterations to use as burnin.
#'
#' @param n.iter the number of interations to run following burnin for the initial jags fit. Defaults to 500 + burnin,
#'
#' @param n.iter.burnin the number of interations to run overall. This should be large. Defaults to 10000.
#'
#' @param n.tries The number of tries to attempt to fit the model and attain good chain mixing. See details below.
#'
#' @param over.disp. If an overdispersed model should be used. Only changes the model fit for poisson and binomial y.type
#' data. For poisson, a negative binomial model will be fit. For binomial a beta model will be fit.
#'
#' @param model The type of model to be fit. Currently takes values of "NEC3param",
#' "NEC4param", "NECsigmoidal", "NECHormesis", "ECx4param", "ECxWeibull1", or "ECxWeibull2".
#'
#' @param init.value.warning Indicates if a warning should be given if the init.fun generated fails.
#'
#' @param sig.val Probability value to use as the lower quantile to test significance of the predictor posterior values
#' against the control, to estimate NEC as an interpolated NOEC value from smooth ECx curves.
#'
#' @details
#'
#' As some concentration-response data will use zero concentration,
#' and there is no distribution on the continuous scale from 0 to in (ie tweedie) available in jags, a small offset
#' is added (1/10^3 of the next lowest value) to zero values of concentration where x.var are gamma distributed.
#'
#' If the initial model fit fails because it cannot be fit by jags (miss-specified priors, invalid
#' initial values), or because there is poor chain mixing fit.jagsNEC will try n.tries many times using the init.fun
#' defined by the write model function. If all those attempts fail, fit.jagsNEC will then try using default values
#' provided by jags. The function will only return an error if all n.tries fail.
#'
#' All models other than "NEC3param" (which is that defined by Fox 2010) are currently undergoing beta testing and are
#' experimental. These should not yet be used for NEC reporting for official purposes. Comments and feedback are welcome,
#' especially reproducible examples of issues, as well as example test cases.
#'
#' All models provide an estimate for NEC. For model types with "NEC" as a prefix, NEC is directly estimated as a paremeter
#' in the model. Models with "ECx" as a prefix are continuos curve models, tyipically used for extracting ECx values
#' from concentration response data. In this instance the NEC value is defined as the concentration at which there is
#' 90 percent certainty (based on the Bayesian posterior estimate) that the response falls below the estimated value of
#' the upper assymptote (top) of the response (i.e the response value is significantly lower than that expected in the case of
#' no exposure).
#'
#' @export
#' @return The $BUGSoutput element of fitted jags model, including an estimate of the NEC value.
#' A posterior sample of the NEC is also available under $sims.list.

fit.jagsNEC <- function(data,
                        x.var,
                        y.var,
                        trials.var = NA,
                        x.type = NA,
                        y.type = NA,
                        burnin = 10000,
                        n.iter = burnin + 500,
                        n.iter.update = 10000,
                        n.tries = 3,
                        params = c("top", "beta", "NEC", "SS", "SSsim"),
                        over.disp = FALSE,
                        model = "NEC3param",
                        init.value.warning = FALSE,
                        added.model = FALSE,
                        sig.val = 0.01,
                        ...) {
  if (added.model == FALSE) {
    data.check <- jagsNEC_input(
      data = data,
      x.var = x.var,
      y.var = y.var,
      trials.var = trials.var,
      x.type = x.type,
      y.type = y.type,
      params = params,
      over.disp = over.disp,
      model = model
    )
    mod.dat <- data.check$mod.dat

    y.type <- data.check$y.type
    x.type <- data.check$x.type
    response <- data.check$response
    data <- data.check$data
    x.dat <- data.check$x.dat
    y.dat <- data.check$y.dat
    params <- data.check$params
    init.fun <- data.check$init.fun
  } else {
    response <- data[, y.var]

    if (y.type == "binomial") {
      mod.dat$trials <- data[, trials.var] # number of "trials"
      response <- data[, y.var] / data[, trials.var]
    }

    y.dat <- data[, y.var]
    x.dat <- data[, x.var]

    mod.file <- Write.jagsModFile(x.type, y.type, mod.dat, params, model)
    init.fun <- mod.file$init.fun
    params <- mod.file$params
  }


  link <- attributes(init.fun)$link

  fit <- fit.jagsMod(mod.dat, init.fun, params, burnin, n.iter, init.value.warning, n.tries)
  J1 <- fit$J1
  all.Js <- fit$all.Js

  if (class(J1) == "try-error") {
    if (length(all.Js) > 0) {
      J1 <- all.Js[[which.min(unlist(lapply(all.Js, FUN = function(x) {
        check.mixing(x)$cv.ratio
      })))]]
      warning(paste(
        "The function fit.jagsNEC was unable to find a set of starting parameters for the ", model,
        " model resulting in good chain mixing. Examine the model using check.chains and carefully evaluate the outcome 
              of the estimated values. Caution should be used in interpreting the model results."
      ), sep = "")
    }
    if (length(all.Js) == 0 & nchar(round(diff(range(x.dat)))) >= 3) {
      stop(paste("The function fit.jagsNEC was unable to fit this model. Your x.var variable ",
        x.var, " covers more than ", nchar(round(diff(range(x.dat)))), " orders of magnitude. 
                 Try fitting the model with x.var on a log scale.",
        se = ""
      ))
    }
    if (length(all.Js) == 0 & nchar(round(diff(range(x.dat)))) < 3) {
      stop("The function fit.jagsNEC was unable to fit this model. Please help us make 
           this software better by emailing a self contained reproducible example to the 
           developers")
    }
  }
  J2 <- update(J1, n.iter = n.iter.update, n.thin = floor((n.iter.update * 0.01)))
  out <- c(J2$BUGSoutput, list(mod.dat = mod.dat, y.type = y.type, x.type = x.type, model = model, link = link))

  min.x <- min(mod.dat$x)
  max.x <- max(mod.dat$x)
  x.seq <- seq(min.x, max.x, length = 1000)

  # extract the relevant model parameters
  extract.params <- c("top", "beta", "NEC", "alpha", "bot", "d", "slope", "EC50")
  extracted.params <- lapply(extract.params, FUN = function(x) {
    quantile(out$sims.list[[x]], c(0.025, 0.5, 0.975))
  })
  names(extracted.params) <- extract.params

  top <- extracted.params$top
  beta <- extracted.params$beta
  NEC <- extracted.params$NEC
  alpha <- extracted.params$alpha
  bot <- extracted.params$bot
  d <- extracted.params$d
  slope <- extracted.params$slope
  EC50 <- extracted.params$EC50

  # set values for unused parameters
  if (is.na(alpha[1])) {
    alpha <- rep(0, 3)
  }

  if (is.na(bot[1])) {
    bot <- rep(0, 3)
    names(bot) <- c("2.5%", "50%", "97.5%")
  }

  if (is.na(d[1])) {
    d <- rep(1, 3)
    names(d) <- c("2.5%", "50%", "97.5%")
  }

  if (is.na(slope[1])) {
    slope <- rep(0, 3)
    names(slope) <- c("2.5%", "50%", "97.5%")
  }

  if (is.na(EC50[1])) {
    if (y.type == "gaussian") {
      EC50 <- extract_ECx.jagsNECfit(out, ECx.val = 50, prob.vals = c(0.025, 0.5, 0.975), type = "relative")
    } else {
      EC50 <- extract_ECx.jagsNECfit(out, ECx.val = 50, prob.vals = c(0.025, 0.5, 0.975))
    }
  }

  # calculate the predicted values based on the median parameter estimates

  if (length(grep("ECx", model)) > 0) {
    mod.class <- "ECx"
  } else {
    mod.class <- "NEC"
  }

  if (mod.class != "ECx") {
    y.pred.m <- predict_NECmod(
      x.vec = x.seq,
      NEC = NEC["50%"], top = top["50%"],
      beta = beta["50%"], d = d["50%"],
      bot = bot["50%"], slope = slope["50%"]
    )
    predicted.y <- predict_NECmod(
      x.vec = mod.dat$x,
      NEC = NEC["50%"], top = top["50%"],
      beta = beta["50%"], d = d["50%"],
      bot = bot["50%"], slope = slope["50%"]
    )
  }
  if (model == "ECx4param") {
    y.pred.m <- predict_ECxmod(
      x.vec = x.seq,
      EC50 = EC50["50%"], top = top["50%"], beta = beta["50%"],
      bot = bot["50%"]
    )
    predicted.y <- predict_ECxmod(
      x.vec = mod.dat$x,
      EC50 = EC50["50%"], top = top["50%"], beta = beta["50%"],
      bot = bot["50%"]
    )
  }

  if (model == "ECxWeibull1") {
    y.pred.m <- predict_WB1mod(
      x.vec = x.seq,
      EC50 = EC50["50%"], top = top["50%"], beta = beta["50%"],
      bot = bot["50%"]
    )
    predicted.y <- predict_WB1mod(
      x.vec = mod.dat$x,
      EC50 = EC50["50%"], top = top["50%"], beta = beta["50%"],
      bot = bot["50%"]
    )
  }

  if (model == "ECxWeibull2") {
    y.pred.m <- predict_WB1mod(
      x.vec = x.seq,
      EC50 = EC50["50%"], top = top["50%"], beta = beta["50%"],
      bot = bot["50%"]
    )
    predicted.y <- predict_WB1mod(
      x.vec = mod.dat$x,
      EC50 = EC50["50%"], top = top["50%"], beta = beta["50%"],
      bot = bot["50%"]
    )
  }

  if (model == "ECxLinear") {
    y.pred.m <- predict_Linearmod(x.vec = x.seq, top = top["50%"], beta = beta["50%"], link = link)
    predicted.y <- predict_Linearmod(x.vec = mod.dat$x, top = top["50%"], beta = beta["50%"], link = link)
  }

  if (model == "ECxExp") {
    y.pred.m <- predict_ECxExpmod(x.vec = x.seq, top = top["50%"], beta = beta["50%"])
    predicted.y <- predict_ECxExpmod(x.vec = mod.dat$x, top = top["50%"], beta = beta["50%"])
  }

  if (model == "ECxsigmoidal") {
    y.pred.m <- predict_ECxsigmoidalmod(x.vec = x.seq, top = top["50%"], beta = beta["50%"], d = d["50%"])
    predicted.y <- predict_ECxsigmoidalmod(x.vec = mod.dat$x, top = top["50%"], beta = beta["50%"], d = d["50%"])
  }

  # calculate the predicted values using the entire posterior
  pred.vals <- c(predict_NECbugsmod(X = out, precision = 1000), list(y.m = y.pred.m))

  # calculate the residuals
  residuals <- response - predicted.y

  # Extract the overdispersion estimate
  od <- mean(out$sims.list$SS > out$sims.list$SSsim)

  # calculate the NEC from the predicted values for the ECx model
  if (mod.class == "ECx") {
    reference <- quantile(pred.vals$posterior[1, ], sig.val)
    out$sims.list$NEC <- sapply(1:out$n.sims, function(x, pred.vals, reference) {
      pred.vals$x[which.min(abs(pred.vals$posterior[, x] - reference))]
    },
    pred.vals = pred.vals, reference = reference
    )

    NEC <- quantile(out$sims.list$NEC, c(0.025, 0.5, 0.975))
  }

  # Put everyting in a list for output
  if (class(out) != "try-error") {
    out <- c(out, list(
      pred.vals = pred.vals,
      NEC = NEC,
      top = top,
      beta = beta,
      alpha = alpha,
      bot = bot,
      d = d,
      EC50 = EC50,
      params = params,
      over.disp = od,
      all.Js = all.Js,
      predicted.y = predicted.y,
      residuals = residuals,
      link = link
    ))

    # assign a class to the output
    class(out) <- "jagsNECfit"
  }

  message(paste("Response variable ", y.var, " modelled using a ", y.type, " distribution.", sep = ""))
  return(out)
}
