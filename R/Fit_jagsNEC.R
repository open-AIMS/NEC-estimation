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
#' @param trials.var the column heading indicating the column for the number of "trials" for binomial response data. If not supplied, the model may run but will not be the model you intended!
#' 
#' @param x.type the statistical distribution to use for the x (concentration) data. This will be guess based on the characteristic of the input data if not supplied. As some concentration-response data will use zero concentration, and there is no distribution on the continuous scale from 0 to in (ie tweedie) available in jags, a small offset is added (1/10^3 of the next lowest value) to zero values of concentration where these are gamma distributed.
#'
#' @param y.type the statistical distribution to use for the y (response) data. This may currently be one of  'binomial', 'poisson',' 'gaussian', or 'gamma'. Others can be added as required, please contact the package maintainer. If not supplied, the appropriate distribution will be guessed based on the distribution of the input data.
#'
#' @param  params A vector of names indicating the parameters that to trace during the jags fit. For the NEC jags model this is typically 'NEC','top' and 'beta'. If left NA, fit.jagsNEC will supply this based on the selected y.type and x.type.
#'
#' @param burnin the number of iterations to use as burnin. 
#' 
#' @param n.iter the number of interations to run following burnin for the initial jags fit. Defaults to 500 + burnin, but small values are good when the model may have trouble fitting and is updated by n.iter.update iterations anyway.
#'
#' @param n.iter.burnin the number of interations to run overall. This should be large. Defaults to 10000.
#'
#' @param n.tries. If the initial model fit fails because it cannot be fit by jags (miss-specified priors, invalid initial values), or because there is poor chain mixing fit.jagsNEC will try n.tries many times using the init.fun defined by write.jags.NECmod. If all those attempts fail, fit.jagsNEC will then try using default values provided by jags. The function will only return an error if all n.tries fail.
#'
#' @export
#' @return The $BUGSoutput element of fitted jags model.

fit.jagsNEC <- function(data,
                        x.var,
                        y.var,
                        trials.var = NA,
                        x.type = NA, 
                        y.type = NA,
                        params = NA,
                        burnin = 5000,
                        n.iter = burnin+500,
                        n.iter.update = 10000,
                        n.tries=10,
                        ...){
  

  
  
  
  # check variable type x.var
  if(is.na(x.type)==T){ # if the x.var is not specified, then guess
    x.dat <- data[,x.var]
    if(class(x.dat)=="integer"){
     stop("jagsNEC does not currently support integer concentration data. Please provide
         a numeric x.var")}    
    if(class(x.dat)=="numeric" & max(x.dat)>1 & min(x.dat)>=0){x.type="gamma"}
    if(class(x.dat)=="numeric" & max(x.dat)<=1 & min(x.dat)>=0){x.type="beta"} 
    if(class(x.dat)=="numeric" & max(x.dat)>1 & min(x.dat)<0){x.type="gaussian"}    
  }

  # check variable type y.var
  if(is.na(y.type)==T){ # if the y.var is not specified, then guess
    y.dat <- data[,y.var]
    if(class(y.dat)=="numeric" & max(y.dat)>1 & min(y.dat)>=0){y.type="gamma"}
    if(class(y.dat)=="numeric" & max(y.dat)<=1 & min(y.dat)>=0){y.type="gamma"} #y.type="beta"
    if(class(y.dat)=="numeric" & min(y.dat)<0){y.type="gaussian"}  
    if(class(y.dat)=="integer" & min(y.dat)>=0 & is.na(trials.var) == TRUE){
      y.type="poisson"} 
    if(is.na(trials.var)!= TRUE & class(y.dat)!="integer"){
      stop("You have supplied a trials.var argument, suggesting you wish to model a binomial.
           Please ensure y.var is an integer representing the number of successes.")}
      
    if(class(y.dat)=="integer" & min(y.dat)>=0 & is.na(trials.var)!= TRUE){
      y.type="binomial"}   
  }   
  
  params <- c("top", "beta", "NEC")
  
  if(y.type=="gamma"){params=c(params,"shape")}
  if(y.type=="gaussian"){params=c(params,"alpha","sigma")}
    
  # error catching for 0 for gamma by adding very small value (no tweedie available in jags)
  if(min(data[,x.var])==0 & x.type=="gamma"){
   tt <- data[,x.var]
   min.val <- min(tt[which(tt>0)])
   data[which(tt==0),x.var] <- tt[which(tt==0)]+(min.val/10^2) 
  } 
  
  if(min(data[,y.var])==0 & y.type=="gamma"){
    tt <- data[,y.var]
    min.val <- min(tt[which(tt>0)])
    data[which(tt==0),y.var] <- tt[which(tt==0)]+(min.val/10^2) 
  } 
  # error catching for 0 for beta by adding very small value (beta does not take zero)
  if(min(data[,x.var])==0 & x.type=="beta"){
    tt <- data[,x.var]
    min.val <- min(tt[which(tt>0)])
    data[which(tt==0),x.var] <- tt[which(tt==0)]+(min.val/10^2) 
  } 
  
  if(min(data[,y.var])==0 & y.type=="beta"){
    tt <- data[,y.var]
    min.val <- min(tt[which(tt>0)])
    data[which(tt==0),y.var] <- tt[which(tt==0)]+(min.val/10^2) 
  } 
  
  # error catching for 1 for beta by subtracting very small value (beta does not take 1)
  if(min(data[,x.var])==1 & x.type=="beta"){
    tt <- data[,x.var]
    max.val <- max(tt[which(tt<1)])
    data[which(tt==1),x.var] <- tt[which(tt==1)]-((max.val)/10^2) 
  } 
  
  if(min(data[,y.var])==1 & y.type=="beta"){
    tt <- data[,y.var]
    max.val <- max(tt[which(tt<1)])
    data[which(tt==1),y.var] <- tt[which(tt==1)]-(max.val/10^2) 
  }
  
  # create jags model data list
  mod.dat <<- list(
    x = data[,x.var],   # concentration
    y = data[,y.var], # response (successes)
    N = nrow(data))  # Sample size
 
  if(y.type=="binomial"){
   mod.dat$trials = data[,trials.var] # number of "trials"
  }
   
  init.fun <- write.jags.NECmod(x=x.type,y=y.type, mod.dat=mod.dat)
  
  all.Js <- list()
  
  J1 <- try(jags(data       = mod.dat,
             inits      = init.fun,
             parameters = params,
             model      = "NECmod.txt",
             n.thin     = 10,
             n.chains   = 5,
             n.burnin   = burnin,
             n.iter     = n.iter), silent = T)
  if(class(J1)!="try-error"){ # make sure the fitted model had good mixing
    all.Js <- c(all.Js, list(J1))   
    if(max(unlist(check.mixing(J1)$cv.test))==1){class(J1)="try-error"}
  }
 
  # if the attempt fails try 10 more times
  w <- 1
    while(class(J1)=="try-error" & w <= n.tries){
    w <- w+1      
    J1 <- try(jags(data       = mod.dat,
                   inits      = init.fun,
                   parameters = params,
                   model      = "NECmod.txt",
                   n.thin     = 10,
                   n.chains   = 5,
                   n.burnin   = burnin,
                   n.iter     = n.iter), silent = T)  
    if(class(J1)!="try-error"){ # make sure the fitted model had good mixing
      all.Js <- c(all.Js, list(J1))   
      if(max(unlist(check.mixing(J1)$cv.test))==1){class(J1)="try-error"}
    }
  }
  
  # if the attempt fails try without initial values
  w <- 1
  while(class(J1)=="try-error" & w <= n.tries){
    w <- w+1 
    J1 <- try(jags(data       = mod.dat,
                   parameters = params,
                   model      = "NECmod.txt",
                   n.thin     = 10,
                   n.chains   = 5,
                   n.burnin   = burnin,
                   n.iter     = n.iter), silent = T) 
    if(class(J1)!="try-error"){ # make sure the fitted model had good mixing
      all.Js <- c(all.Js, list(J1))   
      if(max(unlist(check.mixing(J1)$cv.test))==1){class(J1)="try-error"}
    }    
  }
  
  if(class(J1)=="try-error"){
    if(length(all.Js)>0){
      J1 <- all.Js[[ which.min(unlist(lapply(all.Js, FUN=function(x){check.mixing(x)$cv.ratio})))]]
      warning("The function fit.jagsNEC was unable to find a set of starting parameters resulting in good chain mixing. This may be a result of a weakly defined threshold, in which case try calculating an ECx value instead. If the data show a clear NEC pattern that would normally be easy to fit, please send a reproducible example so we can try and improve our function. While a model has been returned, please evaluate the model using check.chains to assess the validity of the fit, as well as criticially evaluate the outcome of the estimated values. Extreme caution should be used in interpreting the model results.")
         } 
    if(length(all.Js)==0){
      stop("The function fit.jagsNEC was unable to fit this model. Please help us make this software better by emailing a self contained reproducible example to the developers")
    } 
  }
  J2  <- update(J1, n.iter = n.iter.update, n.thin = floor((n.iter.update*0.01)))  
  out <- c(J2$BUGSoutput, list(mod.dat=mod.dat, y.type = y.type, x.type = x.type))
    
  NEC <-  quantile(out$sims.list$NEC,c(0.025, 0.5, 0.975))
  top <-  quantile(out$sims.list$top,c(0.025, 0.5, 0.975))
  beta <-  quantile(out$sims.list$beta,c(0.025, 0.5, 0.975)) 
  alpha <- rep(0,3)#rep(0, length(NEC))
  if(y.type=="gaussian"){
         alpha <-  quantile(out$sims.list$lapha,c(0.025, 0.5, 0.975)) 
  }

  min.x <- min(mod.dat$x)
  max.x <- max(mod.dat$x)
  x.seq <- seq(min.x, max.x, length=100)
  
  y.pred.m <- predict_NECmod(x.vec=x.seq, NEC=NEC["50%"], top=top["50%"], beta=beta["50%"]) 
  pred.vals <- c(predict_NECbugsmod(X=out), list(y.m=y.pred.m))  
  
  out <- c(out, list(
     pred.vals = pred.vals,
     NEC = NEC,
     top =top,
     beta = beta,
     alpha = alpha,
     params = params,
     all.Js = all.Js))
  return(out)    
}


