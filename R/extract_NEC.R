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

#' extract_NSEC
#'
#' Extracts the predicted NSEC value as desired from a jagsNSEC or a jagsMANSEC model fit.
#'
#' @param  X a jag model fit as returned by a call to jags from fit.jagsNSEC
#' 
#' @param sig.val the Probability value to use as the lower quantile to test significance of the predictor posterior values
#' against the control, to estimate NSEC as an interpolated NOEC value from smooth ECx curves.
#' 
#' @param precision The number of unique x values over which to find NSEC - large values will make the NSEC estimate more 
#' precise.
#' 
#' @param posterior A logical value indicating if the full posterior sample of calculated NSEC values should be returned 
#' instead of just the median and 95 credible intervals.
#' 
#' @param xform A function to apply to the returned estimated concentration values
#' 
#' @param x.range A range of x values over which to consider extracting NSEC
#' 
#' @param prob.vals A vector indicating the probability values over which to return the estimated NSEC value. Defaults to 0.5 (median) and 0.025 and 0.975 (95 percent credible intervals). 
#' 
#' @export
#' @return A vector containing the estimated NSEC value, including upper and lower 95 percent Credible Interval bounds
#' 
extract_NSEC <- function(X, sig.val=0.01, precision=10000, posterior = FALSE, type="absolute", 
                        hormesis.def = "control", xform=NA, x.range=NA,
                        prob.vals=c(0.5, 0.025, 0.975), link="identity"){
 
 if(class(X)=="jagsNECfit"){
    NSEC <- extract_NSEC.jagsNECfit(X, sig.val=sig.val, precision=precision, 
                                              posterior = posterior, hormesis.def=hormesis.def, 
                                              xform=xform, 
                                              prob.vals=prob.vals)
  }
 if(class(X)== "jagsMANECfit"){
   NSEC <- extract_NSEC.jagsMANECfit(X, sig.val=sig.val, precision=precision, 
                                   posterior = posterior, xform=xform, x.range=x.range,
                                   prob.vals=prob.vals) 
  }
  
  if(exists("NSEC")==FALSE){
     stop("Failed to estimate NSEC value for the supplied object class. Only jagsNECfit and jagsMANECfit classes are suported.")
  }

  return(NSEC) 

}

#' extract_NSEC.jagsNSEC
#'
#' Extracts the predicted NSEC value as desired from a jagsNSEC model fit object
#'
#' @param  X a jag model fit as returned by a call to jags from fit.jagsNSEC
#' 
#' @param sig.val the Probability value to use as the lower quantile to test significance of the predictor posterior values
#' against the control, to estimate NSEC as an interpolated NOEC value from smooth ECx curves.
#' 
#' @param precision The number of unique x values over which to find NSEC.
#' 
#' @param posterior A logical value indicating if the full posterior sample of calculated NSEC values should be returned 
#' instead of just the median and 95 credible intervals.
#' 
#' @param xform A function to apply to the returned estimated concentration values
#' 
#' @param prob.vals A vector indicating the probability values over which to return the estimated NSEC value. 
#' 
#' @export
#' @return A vector containing the estimated NSEC value, including upper and lower 95 percent Credible Interval bounds

extract_NSEC.jagsNECfit <- function(X, sig.val=0.01, precision=10000, posterior = FALSE, 
                                   hormesis.def = "control", x.range=NA,
                                   xform=NA, prob.vals=c(0.5, 0.025, 0.975)){

      if(sig.val>1 | sig.val<0){
        stop("Supplied sig.val is not in the required range. Please supply a probability value between 0 and 1.")
      }   

    
  if(length(grep("NSEC", X$model))>0){mod.class <- "NSEC"}else{mod.class <- "NSEC"}


    label <- paste("NSEC", sig.val, sep="_")
    
    pred.vals <- predict_NECbugsmod(X, precision=precision, x.range=x.range)
    posterior.sample <- pred.vals$posterior
    x.vec <- pred.vals$'x' 
    
    
    if(X$model=="NSECHormesis" & hormesis.def=="max"){ # remove values prior to the NSEC for the hormesis model
      posterior.sample <- do.call("cbind",
        lapply(1:length(X$sims.list$NEC), FUN=function(x){
         NEC.x <- X$sims.list$NEC[x]
         posterior.x <- posterior.sample[,x]
         posterior.x[which(x.vec<NSEC.x)] <- NA
         return(posterior.x)
        }))
    }
    
    if(X$model=="NSECHormesis"& hormesis.def=="control"){ # remove values greater than the control for the hormesis model
      posterior.sample <- do.call("cbind", 
          lapply(1:ncol(posterior.sample), FUN=function(x){
           control.x <- posterior.sample[1,x]                                    
           posterior.x <- posterior.sample[,x] 
           posterior.x[which(posterior.x>=control.x)] <- NA
           return(posterior.x)
           }))
    }   
    

    
    reference <-  quantile(pred.vals$posterior[1, ], sig.val)

    NSEC.out <- apply(posterior.sample, MARGIN=2, FUN=function(y){
        NSEC.y <- reference
        NSEC.x <- x.vec[which.min(abs(y-NSEC.y))]
        return(NSEC.x)
        
     }) 
  
    
    # calculate the quantile values from the posterior?
    NSEC.estimate <- quantile(unlist(NSEC.out), probs=prob.vals)
    names(NSEC.estimate) <- c(label, paste(label, "lw", sep="_"), paste(label, "up", sep="_"))
    
    # if a transformation is required
    if(class(xform)=="function"){
      NSEC.estimate <- xform(NSEC.estimate)
      NSEC.out <- xform(NSEC.out)
    }   
  
   if(posterior==FALSE){
    return(NSEC.estimate)
  }else{
    return(NSEC.out)}
  
  
}

#' extract_NSEC.jagsMANSEC
#'
#' Extracts the predicted NSEC value as desired from a jagsNSEC model fit object
#'
#' @param  X a fitted jagsMANSEC model object, containing a list of jag model fit as returned by a call to jags from 
#' fit.jagsNSEC
#' 
#' @param sig.val the Probability value to use as the lower quantile to test significance of the predictor posterior values
#' against the control, to estimate NSEC as an interpolated NOEC value from smooth ECx curves.
#' 
#' @param precision The number of unique x values over which to find NSEC.
#' 
#' @param posterior A logical value indicating if the full posterior sample of calculated NSEC values 
#' should be returned instead of just the median and 95 credible intervals.
#' 
#' @param xform A function to apply to the returned estimated concentration values
#' 
#' @param prob.vals A vector indicating the probability values over which to return the estimated NSEC value. 
#' 
#' @export
#' @return A vector containing the estimated NSEC value, including upper and lower 95 percent Credible Interval bounds

extract_NSEC.jagsMANECfit <- function(X, sig.val=0.01, precision=10000, posterior = FALSE,
                                     hormesis.def="control", xform=NA, x.range=NA,
                                     prob.vals=c(0.5, 0.025, 0.975)){
  
  NSEC.out <- unlist(sapply(X$success.models, FUN=function(x){
    base::sample(extract_NSEC.jagsNECfit(X$mod.fits[[x]], 
                                        sig.val=sig.val, 
                                        precision=precision, 
                                        posterior = TRUE, 
                                        x.range = x.range), 
                 round(X$n.sims*X$mod.stats[x, "wi"]))
  }))
  
  label <- paste("NSEC", sig.val, sep="_")
  # calculate the quantile values from the posterior
  NSEC.estimate <- quantile(NSEC.out, probs=prob.vals)
  names(NSEC.estimate) <- c(label, paste(label, "lw", sep="_"), paste(label, "up", sep="_"))
  
  # if a transformation is required
  if(class(xform)=="function"){
    NSEC.estimate <- xform(NSEC.estimate)
    NSEC.out <- xform(NSEC.out)
  }   
  
 
  if(posterior==FALSE){
    return(NSEC.estimate)
  }else{
    return(NSEC.out)}
  
  
}
