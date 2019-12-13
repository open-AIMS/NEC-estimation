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
#' @param model The type of model to be fit. Currently takes values of "", Hockey", "4param" or "basic4param". This is in beta and "Hockey" in particular does not work reliably for all test cases.
#'
#' @export
#' @return The $BUGSoutput element of fitted jags model.

fit.jagsNEC <- function(data,
                        x.var,
                        y.var,
                        trials.var = NA,
                        x.type = NA, 
                        y.type = NA,
                        burnin = 5000,
                        n.iter = burnin+500,
                        n.iter.update = 10000,
                        n.tries=10,
                        params=c("top", "beta", "NEC", "SS", "SSsim"),
                        over.disp=FALSE,
                        model="NEC3param",
                        ...){
  
  y.dat <- data[,y.var]
  x.dat <- data[,x.var] 
  
  # check the x data are numeric
  if(class(x.dat)!="numeric"){
    stop(paste("Your indicated x.var column ", x.var," contains data that is class ", class(x.dat),". 
               The function jagsNEC requires the concentration data (argument x.var) to be numeric.",sep=""))    
  } 
  
  # check the data are lower at high x compared to low x (ie the response variable declines within increase in the x)
  if(mean(y.dat[which(x.dat<mean(x.dat))])< mean(y.dat[which(x.dat>mean(x.dat))])){
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
  if(is.na(match(model, c("", "Hockey", "4param", "basic4param")))){
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
  if(max(data[,x.var])==1 & x.type=="beta"){
    tt <- data[,x.var]
    max.val <- max(tt[which(tt<1)])
    data[which(tt==1),x.var] <- tt[which(tt==1)]-((max.val)/10^2) 
  } 
  
  if(max(data[,y.var])==1 & y.type=="beta"){
    tt <- data[,y.var]
    max.val <- max(tt[which(tt<1)])
    data[which(tt==1),y.var] <- tt[which(tt==1)]-(max.val/10^2) 
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
   if(model=="Hockey"){
    init.fun <- write.jags.Hockey.NECmod(x=x.type,y=y.type, mod.dat=mod.dat)
    params <- c(params, "d")
  } 
  if(model=="4param"){
    init.fun <- write.jags.4param.NECmod(x=x.type,y=y.type, mod.dat=mod.dat)
    params <- c(params, "bot")
  }
  
  if(model=="basic4param"){
    init.fun <- write.jags.basic4param.mod(x=x.type,y=y.type, mod.dat=mod.dat)
    params <- setdiff(c(params, "bot", "EC50"), c("NEC", "alpha"))
  }
  
  all.Js <- list()
  
  J1 <- try(jags(data   = mod.dat,
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
  if(class(J1)=="try-error"){ 
   warning("The generated init.fun failed to yield a valid model. Model based on default jags initial values")
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
      warning("The function fit.jagsNEC was unable to find a set of starting parameters 
              resulting in good chain mixing. This may be a result of a weakly defined 
              threshold, in which case try calculating an ECx value instead. If the data 
              show a clear NEC pattern that would normally be easy to fit, please send a 
              reproducible example so we can try and improve our function. 
              While a model has been returned, please evaluate the model using check.chains 
              to assess the validity of the fit, as well as criticially evaluate the outcome 
              of the estimated values. Extreme caution should be used in interpreting the 
              model results.")
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
  x.seq <- seq(min.x, max.x, length=100)
  
  # extract the model paramters - depending on the model types
  NEC <-   rep(min.x,3); names(NEC) <- c("2.5%",  "50%", "97.5%")
  top <-  quantile(out$sims.list$top,c(0.025, 0.5, 0.975))
  beta <-  quantile(out$sims.list$beta,c(0.025, 0.5, 0.975)) 
  alpha <- rep(0,3)#rep(0, length(NEC))
  bot <- rep(0,3); names(bot) <- c("2.5%",  "50%", "97.5%") #rep(0, length(NEC))
  d <- rep(1,3); names(d) <- c("2.5%",  "50%", "97.5%") #rep(0, length(NEC))
  
  if(model!="basic4param"){  
    NEC <- quantile(out$sims.list$NEC,c(0.025, 0.5, 0.975))
  }
  if(y.type=="gaussian" & model=="NEC3param"){
    alpha <-  quantile(out$sims.list$alpha,c(0.025, 0.5, 0.975)) 
  }
  if(y.type=="gaussian" & model=="Hockey"){
    alpha <-  quantile(out$sims.list$alpha,c(0.025, 0.5, 0.975)) 
  }
  if(model=="4param" | model=="basic4param"){
    bot <-  quantile(out$sims.list$bot, c(0.025, 0.5, 0.975)) 
  }
  if(model=="Hockey"){
    d <-  quantile(out$sims.list$d, c(0.025, 0.5, 0.975)) 
  } 

  # calculate the predicted values based on the median parameter estimates
  if(model!="basic4param"){
    y.pred.m <- predict_NECmod(x.vec=x.seq, 
                             NEC=NEC["50%"], top=top["50%"], beta=beta["50%"], d=d["50%"], bot=bot["50%"])
    predicted.y <- predict_NECmod(x.vec=mod.dat$x,
                                  NEC=NEC["50%"], top=top["50%"], beta=beta["50%"], d=d["50%"], bot=bot["50%"])  
  }else{
    EC50 <-  quantile(out$sims.list$EC50, c(0.025, 0.5, 0.975)) 
    y.pred.m <- predict_ECxmod(x.vec=x.seq, 
                             EC50=EC50["50%"], top=top["50%"], beta=beta["50%"], bot=bot["50%"]) 
    predicted.y <- predict_ECxmod(x.vec=mod.dat$x,
                                  EC50=EC50["50%"], top=top["50%"], beta=beta["50%"], bot=bot["50%"])    
  }
  
  # calculate the predicted values using the entire posterior
  pred.vals <- c(predict_NECbugsmod(X=out), list(y.m=y.pred.m))  
  
  # calculate the residuals
  residuals <-  response - predicted.y
  
  # Extract the overdispersion estimate
  od <- mean(out$sims.list$SS > out$sims.list$SSsim)
  
  # Put everyting in a list for output
  out <- c(out, list(
     pred.vals = pred.vals,
     NEC = NEC,
     top = top,
     beta = beta,
     alpha = alpha,
     bot = bot,
     d = d,
     params = params,
     over.disp = od,
     all.Js = all.Js,
     predicted.y = predicted.y,
     residuals = residuals))
  
  # assign a class to the output
  class(out) <- "jagsNECfit"
  
  message(paste("Response variable ", y.var, " modelled using a ", y.type, " distribution.", sep=""))
  return(out)    
}


