  
write.jags.NECmod <- function(x="gamma", y){  
  
   # binomial y; gamma x ----
   if(x=="gamma" & y=="binomial"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top*exp(-beta*(x[i]-NEC)*step((x[i]-NEC)))
        # response is binomial
        y[i]~dbin(theta[i],trials[i])
        }
        
        # specify model priors
        top ~  dunif(0.0001,0.999) 
        beta~dgamma(0.0001,0.0001)
        NEC~dgamma(0.0001,0.0001) #d norm(3, 0.0001) T(0,)
        
        }
        ", fill=TRUE)
    sink()  #Made model in working directory
    
    init.fun <- function(){list(
      top = rbinom(1, round(mean(mod.dat$trials)), quantile(mod.dat$y/mod.dat$trials, probs=0.75))/
        round(mean(mod.dat$trials)), 
      beta = runif(1,0.0001,0.999),
      NEC = rgamma(1,0.2,0.001))}  
   }

    
    # poisson y; gamma x ----
  if(x=="gamma" & y=="poisson"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top*exp(-beta*(x[i]-NEC)*step((x[i]-NEC)))
        # response is poisson
        y[i]~dpois(theta[i])
        }
        
        # specify model priors
        top ~  dgamma(1,0.001) # dnorm(0,0.001) #T(0,) 
        beta~dgamma(0.0001,0.0001)
        NEC~dgamma(0.0001,0.0001) #d norm(3, 0.0001) T(0,)
        
        }
        ", fill=TRUE)
    sink()  #Made model in working directory
    
    init.fun <- function(){list(
      top = rpois(1,max(mod.dat$y)), 
      beta = rgamma(1,0.2,0.001),
      NEC = rgamma(1,0.2,0.001) #
      
    )}
  }
  
  # poisson y; gaussian x----
  if(x=="gaussian" & y=="poisson"){
    sink("NECmod.txt")
    cat("
        model
        {
        # specify model priors
        top ~  dgamma(1,0.001) # dnorm(0,0.001) #T(0,) #
        beta~dgamma(0.0001,0.0001)
        NEC~dnorm(0, 0.0001)
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top*exp(-beta*(x[i]-NEC)*step((x[i]-NEC)))
        # response is poisson
        y[i]~dpois(theta[i])
        }
        }
        ", fill=TRUE)
    sink()  #Made model in working directory
    
    init.fun <- function(){list(
      top = rpois(1,max(mod.dat$y)), 
      beta = rgamma(1,0.2,0.001),
      NEC = rnorm(1, 0, 2))}
  }
    
  # poisson y; beta x----
  if(x=="beta" & y=="poisson"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top*exp(-beta*(x[i]-NEC)*step((x[i]-NEC)))
        # response is poisson
        y[i]~dpois(theta[i])
        }
        
        # specify model priors
        top ~  dgamma(1,0.001) # dnorm(0,0.001) #T(0,) #
        beta ~ dgamma(0.0001,0.0001)
        NEC ~  dunif(0.0001,0.9999) #dbeta(1,1)
        }
        ", fill=TRUE)
    sink()  #Made model in working directory
    
    init.fun <- function(){list(
      top = rpois(1,max(mod.dat$y)), #rnorm(1,0,1),#
      beta = rgamma(1,0.2,0.001),
      NEC = runif(1, 0.01, 0.9))}
  }
    
  # gamma y; beta x -----
  if(x=="beta" & y=="gamma"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top*exp(-beta*((x[i])-NEC)*step(((x[i])-NEC)))
        # response is gamma
        y[i]~dgamma(shape, shape / (theta[i]))
        }
        
        # specify model priors
        top ~  dgamma(1,0.001) # dnorm(0,0.001) #T(0,) #
        beta ~ dgamma(0.0001,0.0001)
        NEC ~  dunif(0.001,0.999) #dbeta(1,1)
        shape ~ dunif(0,1000)
        }
        ", fill=TRUE)
    sink()  #Made model in working directory
    
    init.fun <- function(){list(
      top = rlnorm(1,log(quantile(mod.dat$y,probs = 0.75)),0.1),
      beta = rgamma(1,0.2,0.001),
      shape = runif(1,0,10),
      NEC = runif(1, 0.3, 0.6))}
  }
    
  # gamma y; gaussian x -----
  if(x=="gaussian" & y=="gamma"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top*exp(-beta*((x[i])-NEC)*step(((x[i])-NEC)))
        # response is gamma
        y[i]~dgamma(shape, shape / (theta[i]))
        }
        
        # specify model priors
        top ~  dgamma(1,0.001) # dnorm(0,0.001) #T(0,) #
        beta ~ dgamma(0.0001,0.0001)
        NEC ~ dnorm(0, 0.0001)
        shape ~ dunif(0,1000)
        }
        ", fill=TRUE)
    sink()  #Made model in working directory
    
    init.fun <- function(){list(
      top = rlnorm(1,log(quantile(mod.dat$y,probs = 0.75)),0.1),
      beta = rgamma(1,0.2,0.001),
      shape = runif(1,0,10),
      NEC = rnorm(1, 0, 1))}
  }
  return(init.fun)
  
}  
  
