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

#' check.chains
#'
#' Generates a plot of MCMC chains and ACF function for a Jags model output.
#'
#' @param  X a jag model fit as returned by a call to jags from fit.jagsNEC
#'
#' @param name an optional character string indicating the label to be placed at the top of the plotting window
#'
#' @export
#' @return A plot of MCMC chains and ACF diagrams for each element in params.

check.chains <- function(X, name = "", pdf.file = "") {
  if (class(X) == "jagsNECfit") {
    params <- X$params
    x <- X$sims.array
    num.chains <- ncol(x[, , params[1]])

    par(mfrow = c(length(params), 2), mar = c(0, 5, 0.5, 0.5), oma = c(4, 0, 2, 0))
    for (i in 1:length(params)) {
      x1 <- as.vector(x[, , params[i]])
      chain.id <- rep(1:num.chains, each = nrow(x[, , params[i]]))
      num.lags <- length(acf(x[, , params[i]][, 1], plot = FALSE)$lag)

      # plot the chains
      plot(1:nrow(x[, , i]), rep(NA, nrow(x[, , params[i]])),
        xaxt = "n",
        ylim = range(x1), main = "", xlab = "", ylab = params[i]
      )
      if (i == length(params)) {
        axis(side = 1)
      }
      for (k in 1:num.chains) {
        lines(1:nrow(x[, , params[i]]), x1[which(chain.id == k)], col = k)
      }

      # plot the acf
      plot(1:num.lags, rep(NA, num.lags),
        xaxt = "n",
        ylim = c(0, 1), xlab = "lag", ylab = "correlation", main = ""
      )
      if (i == length(params)) {
        axis(side = 1)
      }
      for (j in 1:num.chains) {
        acf.j <- acf(x[, , params[i]][, j], plot = FALSE)
        lines(acf.j$lag, acf.j$acf, col = j)
      }
    }
    mtext(name, side = 3, outer = T)
  }
  if (class(X) == "jagsMANECfit") {
    if (nchar(pdf.file) > 0) {
      pdf(file = paste(pdf.file, ".pdf", sep = ""), onefile = TRUE)
    }

    for (m in 1:length(X$mod.fits)) {
      X.m <- X$mod.fits[[m]]
      params <- X.m$params
      x <- X.m$sims.array
      num.chains <- ncol(x[, , params[1]])
      if (nchar(pdf.file) == 0) {
        #x11()
      }
      par(mfrow = c(length(params), 2), mar = c(0, 5, 0.5, 0.5), oma = c(4, 0, 2, 0))
      for (i in 1:length(params)) {
        x1 <- as.vector(x[, , params[i]])
        chain.id <- rep(1:num.chains, each = nrow(x[, , params[i]]))
        num.lags <- length(acf(x[, , params[i]][, 1], plot = FALSE)$lag)

        # plot the chains
        plot(1:nrow(x[, , i]), rep(NA, nrow(x[, , params[i]])),
          xaxt = "n",
          ylim = range(x1), main = "", xlab = "", ylab = params[i]
        )
        if (i == length(params)) {
          axis(side = 1)
        }
        for (k in 1:num.chains) {
          lines(1:nrow(x[, , params[i]]), x1[which(chain.id == k)], col = k)
        }

        # plot the acf
        plot(1:num.lags, rep(NA, num.lags),
          xaxt = "n",
          ylim = c(0, 1), xlab = "lag", ylab = "correlation", main = ""
        )
        if (i == length(params)) {
          axis(side = 1)
        }
        for (j in 1:num.chains) {
          acf.j <- acf(x[, , params[i]][, j], plot = FALSE)
          lines(acf.j$lag, acf.j$acf, col = j)
        }
      }
      mtext(names(X$mod.fits)[m], side = 3, outer = T)
    }
    if (nchar(pdf.file) > 0) {
      dev.off()
    }
  }
}
