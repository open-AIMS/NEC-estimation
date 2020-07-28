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



#' predict_NECmod
#'
#' Calculates predicted y (response) values for a supplied vector of x (concentration) values
#'
#' @param  x.vec the x vector over which to calculate
#'
#' @param NEC the No-Effect-Concentration
#'
#' @param top the estimated value of y (response) over the range of x for which no effect occurs
#'
#' @param beta the exponential decay rate
#'
#' @param alpha the offset of a gaussian y response variable
#'
#' @param d the exponent used in the sigmoidal models
#'
#' @param bot the lower plateau
#'
#' @param slope the slope of the positive linear increase in the hormesis model
#'
#' @export
#' @return A list containing x and fitted y, with up and lw values

predict_NECmod <- function(x.vec, NEC = min(x.vec), top, beta, alpha = 0, d = 1, bot = 0, slope = 0) {
  pre.index <- which(x.vec <= NEC)
  post.index <- which(x.vec > NEC)
  x.seq.pre <- x.vec[pre.index]
  x.seq.post <- x.vec[post.index]

  y.pred <- rep(NA, length(x.vec))

  y.pred.pre <- (top + slope * x.seq.pre) - alpha
  y.pred.post <- bot + (((top + slope * x.seq.post) - bot) * exp(-beta * (x.seq.post - NEC)^d)) - alpha

  y.pred[pre.index] <- y.pred.pre
  y.pred[post.index] <- y.pred.post

  return(y.pred)
}

#' predict_ECxmod
#'
#' Calculates predicted y (response) values for a supplied vector of x (concentration) values
#'
#' @param  x.vec the x vector over which to calculate
#'
#' @param EC50 the 50 percent effect concentration
#'
#' @param top the upper plateau
#'
#' @param beta the exponential decay rate (hillslope)
#'
#' @param bot the lower plateau
#'
#' @export
#' @return A list containing x and fitted y, with up and lw values

predict_ECxmod <- function(x.vec, EC50, top, beta, bot = 0) {
  x.seq <- x.vec
  y.pred <- top + (bot - top) / (1 + exp((EC50 - x.seq) * beta))

  return(y.pred)
}

#' predict_WB1mod
#'
#' Calculates predicted y (response) values for a supplied vector of x (concentration) values
#'
#' @param  x.vec the x vector over which to calculate
#'
#' @param EC50 the 50 percent effect concentration
#'
#' @param top the upper plateau
#'
#' @param beta the exponential decay rate (hillslope)
#'
#' @param bot the lower plateau
#'
#' @export
#' @return A list containing x and fitted y, with up and lw values

predict_WB1mod <- function(x.vec, EC50, top, beta, bot = 0) {
  x.seq <- x.vec
  y.pred <- bot + (top - bot) * exp(-exp(beta * (x.seq - EC50))) # WB1

  return(y.pred)
}

#' predict_WB2mod
#'
#' Calculates predicted y (response) values for a supplied vector of x (concentration) values
#'
#' @param  x.vec the x vector over which to calculate
#'
#' @param EC50 the 50 percent effect concentration
#'
#' @param top the upper plateau
#'
#' @param beta the exponential decay rate (hillslope)
#'
#' @param bot the lower plateau
#'
#' @export
#' @return A list containing x and fitted y, with up and lw values

predict_WB2mod <- function(x.vec, EC50, top, beta, bot = 0) {
  x.seq <- x.vec
  y.pred <- bot + (top - bot) * (1 - exp(-exp(beta * (x.seq - EC50)))) # WB2

  return(y.pred)
}


#' predict_Linearmod
#'
#' Calculates predicted y (response) values for a supplied vector of x (concentration) values
#'
#' @param  x.vec the x vector over which to calculate
#'
#' @param top the upper plateau
#'
#' @param beta the linear decay rate
#'
#' @param link the link function used
#'
#' @export
#' @return A list containing x and fitted y, with up and lw values

predict_Linearmod <- function(x.vec, top, beta, link) {
  x.seq <- x.vec
  y.pred <- xform_link(top - beta * x.seq, use.link = link)

  return(y.pred)
}

#' predict_EcxExpmod
#'
#' Calculates predicted y (response) values for a supplied vector of x (concentration) values
#'
#' @param  x.vec the x vector over which to calculate
#'
#' @param top the upper plateau
#'
#' @param beta the exponential decay rate (hillslope)
#'
#' @export
#' @return A list containing x and fitted y, with up and lw values

predict_ECxExpmod <- function(x.vec, top, beta) {
  x.seq <- x.vec
  y.pred <- top * exp(-beta * (x.seq))

  return(y.pred)
}

#' predict_Ecxsigmoidalmod
#'
#' Calculates predicted y (response) values for a supplied vector of x (concentration) values
#'
#' @param  x.vec the x vector over which to calculate
#'
#' @param top the upper plateau
#'
#' @param beta the exponential decay rate (hillslope)
#'
#' @export
#' @return A list containing x and fitted y, with up and lw values

predict_ECxsigmoidalmod <- function(x.vec, top, beta, d) {
  x.seq <- x.vec
  y.pred <- top * exp(-beta * (x.seq)^d)

  return(y.pred)
}

#' xform_link
#'
#' Calculates predicted y (response) values for a supplied vector of x (concentration) values
#'
#' @param  x.vec the x vector over which to calculate
#'
#' @param link the link function used
#'
#' @export
#' @return A list containing x and fitted y, with up and lw values

xform_link <- function(pred.vec, use.link) {
  y.pred <- pred.vec

  if (use.link == "identity") {
    y.pred <- pred.vec
  }
  if (use.link == "logit") {
    y.pred <- exp(pred.vec) / (1 + exp(pred.vec))
  }
  if (use.link == "log") {
    y.pred <- exp(pred.vec)
  }

  return(y.pred)
}


#' predict_NECbugsmod
#'
#'
#' @param X the jagsNEC model fit (as returned by fit.jagsNEC)
#'
#' @param precision the number of x values over which to predict values
#'
#' @param x.range The range of x values over which to make predictions
#'
#' @export
#' @return A list containing x and fitted y, with up and lw values

predict_NECbugsmod <- function(X, precision = 100, x.range = NA) {
  mod.dat <- X$mod.dat
  min.x <- min(mod.dat$x)
  max.x <- max(mod.dat$x)

  if (is.na(x.range[1])) {
    x.seq <- seq(min.x, max.x, length = precision)
  } else {
    x.seq <- seq(min(x.range), max(x.range), length = precision)
  }


  # for the original model
  if (X$y.type != "gaussian" & X$model == "NEC3param") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_NECmod(
        x.vec = x.seq,
        NEC = X$sims.list$NEC[x],
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x]
      )
    }))
  }

  if (X$y.type == "gaussian" & X$model == "NEC3param") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_NECmod(
        x.vec = x.seq,
        NEC = X$sims.list$NEC[x],
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x],
        alpha = X$sims.list$alpha[x]
      )
    }))
  }

  # for the NECsigmoidal model
  if (X$y.type != "gaussian" & X$model == "NECsigmoidal") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_NECmod(
        x.vec = x.seq,
        NEC = X$sims.list$NEC[x],
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x],
        d = X$sims.list$d[x]
      )
    }))
  }

  if (X$y.type == "gaussian" & X$model == "NECsigmoidal") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_NECmod(
        x.vec = x.seq,
        NEC = X$sims.list$NEC[x],
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x],
        alpha = X$sims.list$alpha[x],
        d = X$sims.list$d[x]
      )
    }))
  }

  # for the 4param model
  if (X$model == "NEC4param") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_NECmod(
        x.vec = x.seq,
        NEC = X$sims.list$NEC[x],
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x],
        bot = X$sims.list$bot[x]
      )
    }))
  }

  if (X$model == "ECx4param") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_ECxmod(
        x.vec = x.seq,
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x],
        EC50 = X$sims.list$EC50[x],
        bot = X$sims.list$bot[x]
      )
    }))
  }

  if (X$model == "ECxWeibull1") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_WB1mod(
        x.vec = x.seq,
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x],
        EC50 = X$sims.list$EC50[x],
        bot = X$sims.list$bot[x]
      )
    }))
  }

  if (X$model == "ECxWeibull2") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_WB2mod(
        x.vec = x.seq,
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x],
        EC50 = X$sims.list$EC50[x],
        bot = X$sims.list$bot[x]
      )
    }))
  }

  # for the NECHormesis model
  if (X$y.type != "gaussian" & X$model == "NECHormesis") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_NECmod(
        x.vec = x.seq,
        NEC = X$sims.list$NEC[x],
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x],
        slope = X$sims.list$slope[x]
      )
    }))
  }

  if (X$y.type == "gaussian" & X$model == "NECHormesis") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_NECmod(
        x.vec = x.seq,
        NEC = X$sims.list$NEC[x],
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x],
        alpha = X$sims.list$alpha[x],
        slope = X$sims.list$slope[x]
      )
    }))
  }

  if (X$model == "ECxLinear") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_Linearmod(
        x.vec = x.seq,
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x],
        link = X$link
      )
    }))
  }

  if (X$model == "ECxExp") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_ECxExpmod(
        x.vec = x.seq,
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x]
      )
    }))
  }

  if (X$model == "ECxsigmoidal") {
    pred.vals.out <- do.call("cbind", lapply(1:X$n.sims, FUN = function(x) {
      predict_ECxsigmoidalmod(
        x.vec = x.seq,
        top = X$sims.list$top[x],
        beta = X$sims.list$beta[x],
        d = X$sims.list$d[x]
      )
    }))
  }

  m.vals <- apply(pred.vals.out, MARGIN = 1, FUN = quantile, probs = 0.5)
  up.vals <- apply(pred.vals.out, MARGIN = 1, FUN = quantile, probs = 0.975)
  lw.vals <- apply(pred.vals.out, MARGIN = 1, FUN = quantile, probs = 0.025)

  return(list(
    x = x.seq,
    y = m.vals,
    up = up.vals,
    lw = lw.vals,
    posterior = pred.vals.out
  ))
}

#' predict
#'
#' Calculated predicted values for a jagsNEC or a jagsMANEC model fit.
#'
#' @param  X a jagsNEC model fit as returned by a call to jags from fit.jagsNEC or fit.jagsMANEC
#'
#' @param precision The number of unique x values over which to find fitted - large values will make the fitted estimate more
#' precise.
#'
#' @param posterior A logical value indicating if the full posterior sample of calculated fitted values should be returned
#' instead of just the median and 95 credible intervals.
#'
#' @param x.range A range of x values over which to consider extracting fitted
#'
#' @param prob.vals A vector indicating the probability values over which to return the estimated fitted value. Defaults to 0.5 (median) and 0.025 and 0.975 (95 percent credible intervals).
#'
#' @export
#' @return A vector containing the estimated fitted value, including upper and lower 95 percent Credible Interval bounds
#'
predict <- function(X, precision = 100, posterior = FALSE, x.range = NA,
                    prob.vals = c(0.5, 0.025, 0.975), link = "identity") {
  if (class(X) == "jagsNECfit") {
    pred.vals <- predict.jagsNECfit(
      X,
      x.range = x.range,
      precision = precision,
      posterior = posterior,
      prob.vals = prob.vals
    )
  }
  if (class(X) == "jagsMANECfit") {
    pred.vals <- predict.jagsMANECfit(
      X,
      precision = precision,
      posterior = posterior,
      x.range = x.range,
      prob.vals = prob.vals
    )
  }

  if (exists("pred.vals") == FALSE) {
    stop("Failed to estimate predicted values for the supplied object class. Only jagsNECfit and jagsMANECfit classes are suported.")
  }

  return(pred.vals)
}

#' predict.jagsNEC
#'
#' Extracts the predicted fitted value as desired from a jagsNEC model fit obeject
#'
#' @inheritParams predict
#'
#' @export
#' @return A vector containing the estimated fitted value, including upper and lower 95 percent Credible Interval bounds

predict.jagsNECfit <- function(X, precision = 100, posterior = FALSE, x.range = NA, prob.vals = c(0.5, 0.025, 0.975)) {
  fitted.out <- predict_NECbugsmod(X, precision = precision, x.range = x.range)
  posterior.predicted <- fitted.out$posterior

  x.vec <- fitted.out$"x"

  # calculate the quantile values from the posterior
  fitted.summary <- data.frame(t(apply(posterior.predicted, MARGIN = 1, FUN = quantile, probs = prob.vals)))
  colnames(fitted.summary) <- c("fit", "lw", "up")
  fitted.summary$x <- x.vec

  if (posterior == FALSE) {
    return(fitted.summary)
  } else {
    return(fitted.out)
  }
}

#' predict.jagsMANEC
#'
#' Extracts the predicted fitted value as desired from a jagsNEC model fit obeject
#'
#' @inheritParams predict
#'
#' @export
#' @return A vector containing the estimated fitted value, including upper and lower 95 percent Credible Interval bounds

predict.jagsMANECfit <- function(X, precision = 100, posterior = FALSE, x.range = NA, prob.vals = c(0.5, 0.025, 0.975)) {
  mod.dat <- X$mod.dat
  min.x <- min(mod.dat$x)
  max.x <- max(mod.dat$x)

  if (is.na(x.range[1])) {
    x.seq <- seq(min.x, max.x, length = precision)
  } else {
    x.seq <- seq(min(x.range), max(x.range), length = precision)
  }

  posterior.predicted <- data.frame(
    lapply(X$success.models, FUN = function(x) {
      posterior.predicted.m <- predict.jagsNECfit(X$mod.fits[[x]],
        precision = precision,
        posterior = TRUE,
        x.range = x.range
      )$posterior
      return(posterior.predicted.m[, base::sample(
        1:X$n.sims,
        round(X$n.sims * X$mod.stats[x, "wi"])
      )])
    })
  )



  # calculate the quantile values from the posterior
  fitted.summary <- data.frame(t(apply(posterior.predicted, MARGIN = 1, FUN = quantile, probs = prob.vals)))
  colnames(fitted.summary) <- c("fit", "lw", "up")
  fitted.summary$x <- x.seq

  # construct the fitted out to match that contained in predict.jagsNEC
  fitted.out <- list(x = x.seq, y = fitted.summary$fit, up = fitted.summary$up, lw = fitted.summary$lw, posterior = posterior.predicted)


  if (posterior == FALSE) {
    return(fitted.summary)
  } else {
    return(fitted.out)
  }
}
