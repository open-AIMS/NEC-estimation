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
#' @param model A vector of the names of model types to be fit. Currently defaults to all availables, including "NEC3param",  
#' "NEC4param", "NECsigmoidal", "NECHormesis" or "ECx4param". This method is under development and desting and should not yet be used for NEC reporting.
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
                        model.set=c("NEC3param", "NEC4param", "NECsigmoidal", "NECHormesis", "ECx4param"),
                        ...){
  
 # Fit each of the models
 mod.fits <- list()   
 for(m in 1: length(model.set)){
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
      silent=T)
    mod.fits <- c(mod.fits, list(fit.m))
    }

  names(mod.fits) <- model.set
  
   # extract model parameters that do not vary across models
   success.models <- model.set[unlist(lapply(mod.fits, FUN=class))=="jagsNECfit"]
   if(length(success.models)==0){
     stop("None of the models fit successfully, 
     try using fit.jagsNEC instead using the default settings as a starting point for trouble shooting.")}else{
       warning(paste("successfully fitted the models: ", paste(success.models, collapse=" ")))
     }
  
  # extract the model statistics for each fit
  mod.stats <- data.frame(DIC = unlist(lapply(mod.fits, FUN=function(x){x$DIC})))
  mod.stats$DIC.delta <- mod.stats$DIC-min(mod.stats$DIC, na.rm = TRUE)
  mod.stats$wi <- wi(mod.stats$DIC)
  mod.stats$pD <- unlist(lapply(mod.fits, FUN=function(x){x$pD}))
  mod.stats$over.disp <- unlist(lapply(mod.fits, FUN=function(x){x$over.disp}))

  mcmc.stats <- unique(do.call("rbind", lapply(mod.fits, FUN=function(x){
               x[c("n.chains", "n.iter", "n.burnin", "n.thin", "n.keep", "n.sims", "y.type", "x.type")]
    })))
  mcmc.list <- do.call("list", mcmc.stats)
  names(mcmc.list) <- colnames(mcmc.stats)
  
  #tt2 <- (lapply(mod.fits, FUN=function(x){
  #  do.call("rbind", x[c("NEC", "top", "beta", "bot", "slope", "d", "alpha")])
  #}))
  
  # model averaged NEC posterior
  NEC.posterior <- unlist(lapply(success.models, FUN=function(x){
    base::sample(mod.fits[[x]]$sims.list$NEC, round(mcmc.list$n.sims*mod.stats[x, "wi"]))
    }))
  
  # model averaged predicted y
  predicted.y <- rowSums(do.call("cbind", lapply(success.models, FUN=function(x){
    mod.fits[[x]]$predicted.y*mod.stats[x, "wi"]
  })))
  
  # model averaged pred.vals
  x <- mod.fits[[success.models[1]]]$pred.vals$x
  
  y <- rowSums(do.call("cbind", lapply(success.models, FUN=function(x){
                      mod.fits[[x]]$pred.vals$y*mod.stats[x, "wi"]
                    })))
  
  up <- rowSums(do.call("cbind", lapply(success.models, FUN=function(x){
                      mod.fits[[x]]$pred.vals$up*mod.stats[x, "wi"]
                    })))
  
  lw <- rowSums(do.call("cbind", lapply(success.models, FUN=function(x){
                      mod.fits[[x]]$pred.vals$lw*mod.stats[x, "wi"]
                    })))
  
  posterior <- Reduce("+", lapply(success.models, FUN=function(x){
                  mod.fits[[x]]$pred.vals$posterior*mod.stats[x, "wi"]    
                     }))
    
  y.m <- rowSums(do.call("cbind", lapply(success.models, FUN=function(x){
                      mod.fits[[x]]$pred.vals$y.m*mod.stats[x, "wi"]
                    })))

  # collate all the elements
  export.list <- 
    c(mcmc.list,
           list(mod.fits=mod.fits,
              mod.dat=mod.dat,
              mod.stats=mod.stats,
              sims.list=list(NEC=NEC.posterior),
              predicted.y=predicted.y,
              residuals=mod.dat$y-predicted.y,
              pred.vals=list(x=x, y=y, up=up, lw=lw, posterior=posterior, y.m=y.m),
              NEC=NEC)
      )
  # assign a class to the output
  class(out) <- "jagsMANECfit"
  
  # return the collated output
  return(export.list)

}