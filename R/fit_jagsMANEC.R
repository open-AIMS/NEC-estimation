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
#' @param  data a data.frame containing the data to use for the model
#' 
#' @param x.var the column heading indicating the concentration (x) variable
#' 
#' @param y.var the column heading indicating the response (y) variable
#' 
#' @param trials.var the column heading indicating the column for the number of "trials" for binomial response data. If not supplied, the model may run but will not be the model you intended!
#' 
#' @param x.type the statistical distribution to use for the x (concentration) data. This will be guess based on the characteristic of the input data if not supplied. As some concentration-response data will use zero concentration, and there is no distribution on the continuous scale from 0 to in (ie tweedie) available in jags, a small offset is added (1/10^3 of the next lowest value) to zero values of concentration where these are gamma distributed.
#'
#' @param y.type the statistical distribution to use for the y (response) data. This may currently be one of  'binomial', 'poisson',' 'gaussian', or 'gamma'. Others can be added as required, please contact the package maintainer. If not supplied, the appropriate distribution will be guessed based on the distribution of the input data.
#'
#' @param  params A vector of names indicating the parameters that to trace during the jags fit. For the NEC jags model this is typically 'NEC','top' and 'beta'. If left out, fit.jagsNEC will supply this based on the selected y.type and x.type.
#'
#' @param burnin the number of iterations to use as burnin. 
#' 
#' @param n.iter the number of interations to run following burnin for the initial jags fit. Defaults to 500 + burnin, but small values are good when the model may have trouble fitting and is updated by n.iter.update iterations anyway.
#'
#' @param n.iter.burnin the number of interations to run overall. This should be large. Defaults to 10000.
#'
#' @param n.tries. If the initial model fit fails because it cannot be fit by jags (miss-specified priors, invalid initial values), or because there is poor chain mixing fit.jagsNEC will try n.tries many times using the init.fun defined by the write model function. If all those attempts fail, fit.jagsNEC will then try using default values provided by jags. The function will only return an error if all n.tries fail.
#'
#' @param over.disp. If an overdispersed model should be used. Only changes the model fit for poisson and binomial y.type data. For poisson, a negative binomial model will be fit. For binomial a beta model will be fit.
#'
#' @param model.set A vector of the names of model types to be fit. Currently defaults to 
#' all available model types. If "NEC" is supplied, only the NEC models will be fit. If "ECx" is supplied,
#'  only continuous curve models will be fit. 
#' 
#' @param sig.val Probability value to use as the lower quantile to test significance of the predictor posterior values
#' against the control, to estimate NEC as an interpolated NOEC value from smooth ECx curves.
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
                        n.iter = burnin+500,
                        n.iter.update = 10000,
                        n.tries=3,
                        params=c("top", "beta", "NEC", "SS", "SSsim"),
                        over.disp=FALSE,
                        model.set="all",
                        sig.val=0.01,
                        ...){
  
  if(model.set=="NEC"){model.set=c("NEC3param", "NEC4param", "NECHormesis", "NECsigmoidal")}
  if(model.set=="ECx"){model.set=c("ECx4param", "ECxWeibull1", "ECxWeibull2","ECxLinear")}
  if(model.set=="all"){model.set=c("NEC3param", "NEC4param", "NECHormesis", "NECsigmoidal", 
                                   "ECxLinear", "ECxExp", "ECxsigmoidal",
                                   "ECx4param", "ECxWeibull1", "ECxWeibull2")}
  if(model.set=="bot_free"){model.set=c("NEC3param", "NECHormesis", "NECsigmoidal", 
                                   "ECxLinear", "ECxExp", "ECxsigmoidal")}
  
 # Fit each of the models
 mod.fits <- vector(mode = 'list', length = length(model.set))
 names(mod.fits) <- model.set
 
 for(m in seq_along(model.set)){
    model <- model.set[m] 
    fit.m <- try(
      fit.jagsNEC(data=data,
                            x.var=x.var,
                            y.var=y.var,
                            trials.var = trials.var,
                            x.type = x.type, 
                            y.type = y.type,
                            burnin = burnin,
                            n.iter = n.iter,
                            n.iter.update = n.iter.update,
                            n.tries=n.tries,
                            params=params,
                            over.disp=over.disp,
                            model=model), 
      silent = TRUE)
    if (!inherits(fit.m, 'try-error')) {
      mod.fits[[m]] <- fit.m  
    } else {
      mod.fits[[m]] <- NA 
    }
    
    }

   

  # collate all the elements
  export.list <- c(extract_modStats(mod.fits), 
                   list(data=data, x.var=x.var, y.var=y.var, trials.var=trials.var, over.disp=over.disp))
  # assign a class to the output
  class(export.list) <- "jagsMANECfit"
  
  # return the collated output
  return(export.list)

}
