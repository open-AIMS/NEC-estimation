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
#' @param n.tries. The number of tries to attempt to fit the model and attain good chain mixing. See details below.
#'
#' @param over.disp. If an overdispersed model should be used. Only changes the model fit for poisson and binomial y.type 
#' data. For poisson, a negative binomial model will be fit. For binomial a beta model will be fit.
#'
#' @param model The type of model to be fit. Currently takes values of "NEC3param",  
#' "NEC4param", "NECsigmoidal", "NECHormesis", "ECx4param", "ECxWeibull1", or "ECxWeibull2".
#' 
#' @param init.value.warning Indicates if a warning should be given if the init.fun generated fails.
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
                        n.iter = burnin+500,
                        n.iter.update = 10000,
                        n.tries=3,
                        params=c("top", "beta", "NEC", "SS", "SSsim"),
                        over.disp=FALSE,
                        model="NEC3param",
                        init.value.warning=FALSE,
                        ...){
  
  y.dat <- data[, y.var]
  x.dat <- data[, x.var] 
   
  # check the x data are numeric
  if(class(x.dat)!="numeric"){
    stop(paste("Your indicated x.var column ", x.var," contains data that is class ", class(x.dat),". 
               The function jagsNEC requires the concentration data (argument x.var) to be numeric.",sep=""))    
  } 
  
  # check data contains only finite values
  test.x <- mean(x.dat)
  test.y <- mean(y.dat)
  if(is.finite(test.x)!=TRUE){
    stop("Your x.var column contains values that are not finite.")    
  } 
  if(is.finite(test.y)!=TRUE){
    stop("Your y.var column contains values that are not finite.")    
  } 
  
  # check the data are lower at high x compared to low x (ie the response variable declines within increase in the x)
  if(mean(y.dat[which(x.dat<mean(x.dat))])< mean(y.dat[which(x.dat>mean(x.dat))]) & model != "NECHormesis"){
    stop("The mean value of the response for the lower half of the 
         concentration data are lower than that of the upper half of the concentration data. 
         jagsNEC only fits concentration response data where the 
         response declines with increasing values of concentration.")    
  }

  # check variable type x.var
  if(is.na(x.type)==TRUE){ # if the x.var is not specified, then guess
     if(class(x.dat)=="integer"){
     stop("jagsNEC does not currently support integer concentration data. Please provide
         a numeric x.var")}    
    if(class(x.dat)=="numeric" & max(x.dat)>1 & min(x.dat)>=0){x.type="gamma"}
    if(class(x.dat)=="numeric" & max(x.dat)<=1 & min(x.dat)>=0){x.type="beta"} 
    if(class(x.dat)=="numeric" & max(x.dat)>1 & min(x.dat)<0){x.type="gaussian"}    
  }

  # check variable type y.var
  if(is.na(y.type)==T){ # if the y.var is not specified, then guess
    if(class(y.dat)=="numeric" & max(y.dat)>1 & min(y.dat)>=0){y.type="gamma"}
    if(class(y.dat)=="numeric" & max(y.dat)<=1 & min(y.dat)>=0){y.type="beta"}
    if(class(y.dat)=="numeric" & min(y.dat)<0){y.type="gaussian"}  
    if(class(y.dat)=="integer" & min(y.dat)>=0 & is.na(trials.var) == TRUE){
      y.type="poisson"} 
    if(is.na(trials.var)!= TRUE & class(y.dat)!="integer"){
      stop("You have supplied a trials.var argument, suggesting you wish to model a binomial.
           Please ensure y.var is an integer representing the number of successes.")}
      
    if(class(y.dat)=="integer" & min(y.dat)>=0 & is.na(trials.var)!= TRUE){
      y.type="binomial"}   
  } 
  
  # check there is a valid model type
  if(is.na(match(model, c("NEC3param", "NECsigmoidal", "NEC4param", "NECHormesis",
                          "ECx4param", "ECxWeibull1", "ECxWeibull2")
                 ))){
    stop("The model type you have specified does not extist.")    
  }
  
  if(y.type=="poisson" & over.disp==TRUE){
          y.type="negbin"}
  if(y.type=="binomial" & over.disp==TRUE){
          y.type <- "beta"
          data[,y.var] <-  data[,y.var]/data[,trials.var]  
          
          }  

  if(y.type=="gamma"){params=c(params,"shape")}
  if(y.type=="gaussian"){params=c(params,"alpha","sigma")}
  if(y.type=="negbin"){params=c(params,"size")}  
    
  # error catching for 0 for gamma by adding very small value (no tweedie available in jags)
  if(min(data[,x.var])==0 & x.type=="gamma"){
   tt <- data[,x.var]
   min.val <- min(tt[which(tt>0)])
   data[which(tt==0),x.var] <- tt[which(tt==0)]+(min.val/10) 
  } 
  
  if(min(data[,y.var])==0 & y.type=="gamma"){
    tt <- data[,y.var]
    min.val <- min(tt[which(tt>0)])
    data[which(tt==0),y.var] <- tt[which(tt==0)]+(min.val/10) 
  } 
  # error catching for 0 for beta by adding very small value (beta does not take zero)
  if(min(data[,x.var])==0 & x.type=="beta"){
    tt <- data[,x.var]
    min.val <- min(tt[which(tt>0)])
    data[which(tt==0),x.var] <- tt[which(tt==0)]+(min.val/10) 
  } 
  
  if(min(data[,y.var])==0 & y.type=="beta"){
    tt <- data[,y.var]
    min.val <- min(tt[which(tt>0)])
    data[which(tt==0),y.var] <- tt[which(tt==0)]+(min.val/10) 
  } 
  
  # error catching for 1 for beta by subtracting very small value (beta does not take 1)
  if(max(data[,x.var])==1 & x.type=="beta"){
    tt <- data[,x.var]
    max.val <- max(tt[which(tt<1)])
    data[which(tt==1),x.var] <- tt[which(tt==1)]-((max.val)/10) 
  } 
  
  if(max(data[,y.var])==1 & y.type=="beta"){
    tt <- data[,y.var]
    max.val <- max(tt[which(tt<1)])
    data[which(tt==1),y.var] <- tt[which(tt==1)]-(max.val/10) 
  }
  
  # create jags model data list
  mod.dat <<- list(
    x = data[,x.var],   # concentration
    y = data[,y.var], # response (successes)

    N = nrow(data))  # Sample size

  
  response = data[,y.var]
 
  if(y.type=="binomial"){
   mod.dat$trials = data[, trials.var] # number of "trials"
   response = data[, y.var]/data[,trials.var]
  }
  
  # set the type of model to fit
  if(model=="NEC3param"){
     init.fun <- write.jags.NEC3param.mod(x=x.type,y=y.type, mod.dat=mod.dat) 
  }
   if(model=="NECsigmoidal"){
    init.fun <- write.jags.NECsigmoidal.mod(x=x.type,y=y.type, mod.dat=mod.dat)
    params <- c(params, "d")
  } 
  if(model=="NEC4param"){
    init.fun <- write.jags.NEC4param.mod(x=x.type,y=y.type, mod.dat=mod.dat)
    params <- setdiff(c(params, "bot"), c("alpha"))
  }
  
  if(model=="NECHormesis"){
    init.fun <- write.jags.NECHormesis.mod(x=x.type,y=y.type, mod.dat=mod.dat)
    params <- c(params, "slope")
  }
  
  if(model=="ECx4param"){
    init.fun <- write.jags.ECx4param.mod(x=x.type,y=y.type, mod.dat=mod.dat)
    params <- setdiff(c(params, "bot", "EC50"), c("NEC", "alpha"))
  }
  
  if(model=="ECxWeibull1"){
    init.fun <- write.jags.ECxWeibull1.mod(x=x.type,y=y.type, mod.dat=mod.dat)
    params <- setdiff(c(params, "bot", "EC50"), c("NEC", "alpha"))
  }
  
  if(model=="ECxWeibull2"){
    init.fun <- write.jags.ECxWeibull2.mod(x=x.type,y=y.type, mod.dat=mod.dat)
    params <- setdiff(c(params, "bot", "EC50"), c("NEC", "alpha"))
  }
  
  all.Js <- list()
  
  warn = getOption("warn")
  options(warn=-1)
  J1 <- try(jags(data   = mod.dat,
             inits      = init.fun,
             parameters = params,
             model      = "NECmod.txt",
             n.thin     = 10,
             n.chains   = 5,
             n.burnin   = burnin,
             n.iter     = n.iter,
             progress.bar = "none"), silent = T)
  options(warn=warn)
  
  if(class(J1)!="try-error"){ # make sure the fitted model had good mixing
    all.Js <- c(all.Js, list(J1))   
    if(max(unlist(check.mixing(J1)$cv.test))==1){class(J1)="try-error"}
  }
  if(init.value.warning==TRUE){
  warning("The generated init.fun failed to yield a valid model. Model based on default jags initial values")
  }
  #if the attempt fails try 10 more times
  warn = getOption("warn")
  options(warn=-1)
  w <- 1
    while(class(J1)=="try-error" & w <= n.tries){
    w <- w+1      
    J1 <- try(R2jags::jags(data       = mod.dat,
                   inits      = init.fun,
                   parameters = params,
                   model      = "NECmod.txt",
                   n.thin     = 10,
                   n.chains   = 5,
                   n.burnin   = burnin,
                   n.iter     = n.iter,
                   progress.bar = "none"), silent = T)  
    if(class(J1)!="try-error"){ # make sure the fitted model had good mixing
      all.Js <- c(all.Js, list(J1))   
      if(max(unlist(check.mixing(J1)$cv.test))==1){class(J1)="try-error"}
    }
  }
  
  # if the attempt fails try without initial values
  w <- 1
  while(class(J1)=="try-error" & w <= n.tries){
    w <- w+1 
    J1 <- try(R2jags::jags(data       = mod.dat,
                   parameters = params,
                   model      = "NECmod.txt",
                   n.thin     = 10,
                   n.chains   = 5,
                   n.burnin   = burnin,
                   n.iter     = n.iter,
                   progress.bar = "none"), silent = T) 
    if(class(J1)!="try-error"){ # make sure the fitted model had good mixing
      all.Js <- c(all.Js, list(J1))   
      if(max(unlist(check.mixing(J1)$cv.test))==1){class(J1)="try-error"}
    }    
  }
  #options(warn=warn)
  
  if(class(J1)=="try-error"){
    if(length(all.Js)>0){
      J1 <- all.Js[[ which.min(unlist(lapply(all.Js, FUN=function(x){check.mixing(x)$cv.ratio})))]]
      warning(paste("The function fit.jagsNEC was unable to find a set of starting parameters for the ", model,  
              " model resulting in good chain mixing. Examine the model using check.chains and carefully evaluate the outcome 
              of the estimated values. Caution should be used in interpreting the model results."), sep="")
         } 
    if(length(all.Js)==0 & nchar(round(diff(range(x.dat))))>=3){
      stop(paste("The function fit.jagsNEC was unable to fit this model. Your x.var variable ", 
                 x.var," covers more than ", nchar(round(diff(range(x.dat))))," orders of magnitude. 
                 Try fitting the model with x.var on a log scale.", se=""))
    } 
    if(length(all.Js)==0 & nchar(round(diff(range(x.dat))))<3){
      stop("The function fit.jagsNEC was unable to fit this model. Please help us make 
           this software better by emailing a self contained reproducible example to the 
           developers")
    } 
  }
  J2  <- update(J1, n.iter = n.iter.update, n.thin = floor((n.iter.update*0.01)))  
  out <- c(J2$BUGSoutput, list(mod.dat=mod.dat, y.type = y.type, x.type = x.type, model = model))
    
  min.x <- min(mod.dat$x)
  max.x <- max(mod.dat$x)
  x.seq <- seq(min.x, max.x, length=1000)

  # extract the relevant model parameters
  extract.params <- c("top", "beta", "NEC", "alpha", "bot", "d", "slope", "EC50")
  extracted.params <- lapply(extract.params, FUN=function(x){
    quantile(out$sims.list[[x]],c(0.025, 0.5, 0.975))
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
  if(is.na(alpha[1])){
    alpha <- rep(0,3)}
  
  if(is.na(bot[1])){ 
    bot <- rep(0,3); names(bot) <- c("2.5%",  "50%", "97.5%")}
  
  if(is.na(d[1])){  
    d <- rep(1,3); names(d) <- c("2.5%",  "50%", "97.5%")}
  
  if(is.na(slope[1])){ 
    slope <- rep(0,3); names(slope) <- c("2.5%",  "50%", "97.5%")}
  
  if(is.na(EC50[1])){
    if(y.type !="gaussian"){
      EC50 <- extract_ECx.jagsNECfit(out, ECx.val = 50, prob.vals = c(0.025, 0.5, 0.975))
      }
    if(y.type =="gaussian"){
      EC50 <- extract_ECx.jagsNECfit(out, ECx.val = 50, prob.vals = c(0.025, 0.5, 0.975), type="relative")
      }  
  }
  
  # calculate the predicted values based on the median parameter estimates
 
  if(length(grep("ECx", model))>0){mod.class <- "ECx"}else{ mod.class <- "NEC"}
    
  if(mod.class!="ECx"){
    y.pred.m <- predict_NECmod(x.vec=x.seq, 
                             NEC=NEC["50%"], top=top["50%"],  
                             beta=beta["50%"], d=d["50%"],
                             bot=bot["50%"], slope=slope["50%"])
    predicted.y <- predict_NECmod(x.vec=mod.dat$x,
                                  NEC=NEC["50%"], top=top["50%"], 
                                  beta=beta["50%"], d=d["50%"], 
                                  bot=bot["50%"], slope=slope["50%"])}
  if(model=="ECx4param"){
    y.pred.m <- predict_ECxmod(x.vec=x.seq, 
                             EC50=EC50["50%"], top=top["50%"], beta=beta["50%"], 
                             bot=bot["50%"]) 
    predicted.y <- predict_ECxmod(x.vec=mod.dat$x,
                                  EC50=EC50["50%"], top=top["50%"], beta=beta["50%"], 
                                  bot=bot["50%"])    
  }
  
  if(model=="ECxWeibull1"){
    y.pred.m <- predict_WB1mod(x.vec=x.seq, 
                               EC50=EC50["50%"], top=top["50%"], beta=beta["50%"], 
                               bot=bot["50%"]) 
    predicted.y <- predict_WB1mod(x.vec=mod.dat$x,
                                  EC50=EC50["50%"], top=top["50%"], beta=beta["50%"], 
                                  bot=bot["50%"])    
  }
  
  if(model=="ECxWeibull2"){
    y.pred.m <- predict_WB1mod(x.vec=x.seq, 
                               EC50=EC50["50%"], top=top["50%"], beta=beta["50%"], 
                               bot=bot["50%"]) 
    predicted.y <- predict_WB1mod(x.vec=mod.dat$x,
                                  EC50=EC50["50%"], top=top["50%"], beta=beta["50%"], 
                                  bot=bot["50%"])    
  }
  # calculate the predicted values using the entire posterior
  pred.vals <- c(predict_NECbugsmod(X=out, precision = 1000), list(y.m=y.pred.m))  
  
  # calculate the residuals
  residuals <-  response - predicted.y
  
  # Extract the overdispersion estimate
  od <- mean(out$sims.list$SS > out$sims.list$SSsim)
  
  # calculate the NEC from the predicted values for the ECx model
  if(mod.class=="ECx"){
    reference <-  quantile(pred.vals$posterior[1, ], 0.01)
    out$sims.list$NEC  <-  sapply(1:out$n.sims, function (x, pred.vals, reference) {
      pred.vals$x[which.min(abs(pred.vals$posterior[, x] - reference))]}, 
      pred.vals = pred.vals, reference = reference)

    NEC <- quantile(out$sims.list$NEC, c(0.025, 0.5, 0.975))  
  }
  
  # Put everyting in a list for output
  if(class(out)!="try-error"){
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
     residuals = residuals))
  
  # assign a class to the output
  class(out) <- "jagsNECfit"
  }

  message(paste("Response variable ", y.var, " modelled using a ", y.type, " distribution.", sep=""))
  return(out)    
}


