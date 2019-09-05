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
#' @param trials.var the column heading indicating the column for the number of "trials" for binomial response data
#' 
#' @param x.type the statistical distribution to use for the x (concentration) data. The default is "gamma" as most concentration data are gamma (ie continuous on the scale of >0 to inf). As some experiments will use zero concentration, and there is no distribution on the continuous scale from 0 to in (ie tweedie) available in jags, a small offset is added (1/10^5 of the next lowest value) to zero values of concentration
#'
#' @param y.type the statistical distribution to use for the y (response) data. This may currently be one of  'binomial', 'poisson', or 'gamma'. Others can be added as required, please contact the package maintainer.
#'
#' @param  params A vector of names indicating the parameters that were traced during the jags fit. For the NEC jags model this is typically 'NEC','top' and 'beta', whicha are the defaults. 
#'
#' @param burnin the number of iterations to use as burning.
#' 
#' @param n.iter the number of interations to run overall. Defaults to 2 times burnin.
#'
#' @export
#' @return The $BUGSoutput element of fitted jags model.

fit.jagsNEC <- function(data,
                        x.var,
                        y.var,
                        trials.var = NA,
                        x.type = "gamma", 
                        y.type,
                        params = c("top", "beta", "NEC"), # default params for typical NEC model
                        burnin = 1000,
                        n.iter = burnin*2
                        ){
  # error catching for 0 concentration for gamma by adding very small value (no tweedie available in jags)
  if(min(data[,x.var])==0 & x.type=="gamma"){
   tt <- data[,x.var]
   min.val <- min(tt[which(tt>0)])
   data[which(tt==0),x.var] <- tt[which(tt==0)]+(min.val/10^5) 
  } 

  mod.dat <<- list(
    x = data[,x.var],   # concentration
    y = data[,y.var], # response (successes)
    N = nrow(data))  # Sample size
 
  if(y.type=="binomial"){
   mod.dat$trials = data[,trials.var] # number of "trials"
  }
   
  init.fun <- write.jags.NECmod(x=x.type,y=y.type)
  
  J1 <- jags(data       = mod.dat,
             inits      = init.fun,
             parameters = params,
             model      = "NECmod.txt",
             n.thin     = 10,
             n.chains   = 5,
             n.burnin   = burnin,
             n.iter     = n.iter)
  out <- c(J1$BUGSoutput, list(mod.dat=mod.dat, y.type = y.type))
  
  NEC <-  quantile(out$sims.list$NEC,c(0.025, 0.5, 0.975))
  top <-  quantile(out$sims.list$top,c(0.025, 0.5, 0.975))
  beta <-  quantile(out$sims.list$beta,c(0.025, 0.5, 0.975)) 
  
  min.x <- min(mod.dat$x)
  max.x <- max(mod.dat$x)
  x.seq <- seq(min.x, max.x, length=100)
  
  y.pred.m <- predict.NECmod(x.vec=x.seq, NEC=NEC["50%"], top=top["50%"], beta=beta["50%"]) 
  pred.vals <- c(predict.NECbugsmod(X=out), list(y.m=y.pred.m))  
  
  out <- c(out, list(
     pred.vals = pred.vals,
     NEC = NEC,
     top =top,
     beta = beta))
  return(out)
  
}


