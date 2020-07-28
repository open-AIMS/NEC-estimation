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

#' Write.jagsModFile
#'
#' Writes a jags model files and generates the required init.fun and param vector.
#'
#' @param  x.type The type of predictor data to write the model for
#'
#' @param y.type the type of response data to write the model for
#'
#' @param mod.dat a list containing the jags model data
#'
#' @param params a character vector containing the parameters of the model
#'
#' @param model the type of model to fit
#'
#' @export
#' @return In init.fun and updated params vector

Write.jagsModFile <- function(x.type, y.type, mod.dat, params, model) {
  # set the type of model to fit
  if (model == "NEC3param") {
    init.fun <- write.jags.NEC3param.mod(x = x.type, y = y.type, mod.dat = mod.dat)
  }
  if (model == "NECsigmoidal") {
    init.fun <- write.jags.NECsigmoidal.mod(x = x.type, y = y.type, mod.dat = mod.dat)
    params <- c(params, "d")
  }
  if (model == "NEC4param") {
    init.fun <- write.jags.NEC4param.mod(x = x.type, y = y.type, mod.dat = mod.dat)
    params <- setdiff(c(params, "bot"), c("alpha"))
  }

  if (model == "NECHormesis") {
    init.fun <- write.jags.NECHormesis.mod(x = x.type, y = y.type, mod.dat = mod.dat)
    params <- c(params, "slope")
  }

  if (model == "ECx4param") {
    init.fun <- write.jags.ECx4param.mod(x = x.type, y = y.type, mod.dat = mod.dat)
    params <- setdiff(c(params, "bot", "EC50"), c("NEC", "alpha"))
  }

  if (model == "ECxWeibull1") {
    init.fun <- write.jags.ECxWeibull1.mod(x = x.type, y = y.type, mod.dat = mod.dat)
    params <- setdiff(c(params, "bot", "EC50"), c("NEC", "alpha"))
  }

  if (model == "ECxWeibull2") {
    init.fun <- write.jags.ECxWeibull2.mod(x = x.type, y = y.type, mod.dat = mod.dat)
    params <- setdiff(c(params, "bot", "EC50"), c("NEC", "alpha"))
  }

  if (model == "ECxLinear") {
    init.fun <- write.jags.ECxLinear.mod(x = x.type, y = y.type, mod.dat = mod.dat)
    params <- setdiff(params, c("NEC", "alpha"))
  }

  if (model == "ECxExp") {
    init.fun <- write.jags.ECxExp.mod(x = x.type, y = y.type, mod.dat = mod.dat)
    params <- setdiff(params, c("NEC", "alpha"))
  }

  if (model == "ECxsigmoidal") {
    init.fun <- write.jags.ECxsigmoidal.mod(x = x.type, y = y.type, mod.dat = mod.dat)
    params <- setdiff(c(params, "d"), c("NEC", "alpha"))
  }
  return(list(init.fun = init.fun, params = params))
}
