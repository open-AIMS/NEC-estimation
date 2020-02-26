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

#' write.jags.ECxLinear.mod
#'
#' Writes an EC50 model (log logistic) file and generates a function for initial values to pass to jags
#' 
#' @param x the statistical distribution to use for the x (concentration) data. This may currently be one of  'beta', 'gaussian', or 'gamma'. Others can be added as required, please contact the package maintainer.
#' 
#' @param y the statistical distribution to use for the y (response) data. This may currently be one of  'binomial', 'beta', 'poisson', 'gaussian', or 'gamma'. Others can be added as required, please contact the package maintainer.
#'
#' @export
#' @return an init function to pass to jags

write.jags.ECxLinear.mod <- function(x="gamma", y, mod.dat){  

  # binomial y; gamma x ----
   if(x=="gamma" & y=="binomial"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top - x[i]*beta
        # response is binomial
        y[i]~dbin(theta[i],trials[i])
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        top ~  dunif(0.0001,0.99)
         
        # pearson residuals
        for (i in 1:N) {
         ExpY[i] <- theta[i]*trials[i] 
         VarY[i] <- trials[i]*theta[i] * (1 - theta[i])
         E[i] <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
         ysim[i] ~  dbin(theta[i],trials[i])
         Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
         D[i]    <- E[i]^2     #Squared residuals for original data
         Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])

        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = runif(1,0.0001,0.999),#rlnorm(1,0,1), #
      top =runif(1, 0.01, 0.09))}#rgamma(1,0.2,0.001))}  
   }

  # binomial y; gaussian x ----
  if(x=="gaussian" & y=="binomial"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top - x[i]*beta
        # response is binomial
        y[i]~dbin(theta[i],trials[i])
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        top ~  dunif(0.0001,0.99)
         
        # pearson residuals
        for (i in 1:N) {
         ExpY[i] <- theta[i]*trials[i] 
         VarY[i] <- trials[i]*theta[i] * (1 - theta[i])
         E[i] <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
         ysim[i] ~  dbin(theta[i],trials[i])
         Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
         D[i]    <- E[i]^2     #Squared residuals for original data
         Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])
        
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = runif(1,0.0001,0.999),
      top =runif(1, 0.01, 0.09))} 
  }
  
  # binomial y; beta x ----
  if(x=="beta" & y=="binomial"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top - x[i]*beta
        # response is binomial
        y[i]~dbin(theta[i],trials[i])
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        top ~  dunif(0.0001,0.99)
         
        # pearson residuals
        for (i in 1:N) {
        ExpY[i] <- theta[i]*trials[i] 
        VarY[i] <- trials[i]*theta[i] * (1 - theta[i])
        E[i] <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
        ysim[i] ~  dbin(theta[i],trials[i])
        Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
        D[i]    <- E[i]^2     #Squared residuals for original data
        Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])        
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = runif(1,0.0001,0.999),
      top =runif(1, 0.01, 0.09))} 
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
        theta[i]<-top -x[i]*beta
        # response is poisson
        y[i]~dpois(theta[i])
        }
        
        # specify model priors
        beta~dgamma(0.0001,0.0001)
        top ~  dgamma(1,0.001)
        
        # pearson residuals
        for (i in 1:N) {
         ExpY[i] <- theta[i] 
         VarY[i] <- theta[i]
         E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }

        # overdispersion
        for (i in 1:N) {
         ysim[i] ~  dpois(theta[i])
         Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
         D[i]    <- E[i]^2     #Squared residuals for original data
         Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
         }
         SS    <- sum(D[1:N])
         SSsim <- sum(Dsim[1:N])

        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = runif(1,0.0001,0.999), #rlnorm(1,0,1), #rgamma(1,0.2,0.001),
      top =rpois(1,round(mean(mod.dat$y))))} 
  }
  
  # poisson y; gaussian x----
  if(x=="gaussian" & y=="poisson"){
    sink("NECmod.txt")
    cat("
        model
        {

        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top -x[i]*beta
        # response is poisson
        y[i]~dpois(theta[i])
        }

        # specify model priors
        beta~dgamma(0.0001,0.0001)
        top ~  dgamma(1,0.001)

        # pearson residuals
        for (i in 1:N) {
         ExpY[i] <- theta[i] 
         VarY[i] <- theta[i]
         E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }

        # overdispersion
        for (i in 1:N) {
         ysim[i] ~  dpois(theta[i])
         Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
         D[i]    <- E[i]^2     #Squared residuals for original data
         Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
         }
         SS    <- sum(D[1:N])
         SSsim <- sum(Dsim[1:N])
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = rgamma(1,0.2,0.001),
      top =rpois(1,round(mean(mod.dat$y))))} 
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
        theta[i]<-top -x[i]*beta
        # response is poisson
        y[i]~dpois(theta[i])
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        top ~  dgamma(1,0.001)

        # pearson residuals
        for (i in 1:N) {
         ExpY[i] <- theta[i] 
         VarY[i] <- theta[i]
         E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }

        # overdispersion
        for (i in 1:N) {
         ysim[i] ~  dpois(theta[i])
         Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
         D[i]    <- E[i]^2     #Squared residuals for original data
         Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
         }
         SS    <- sum(D[1:N])
         SSsim <- sum(Dsim[1:N])
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = rgamma(1,0.2,0.001),
      top =rpois(1,round(mean(mod.dat$y))))} 
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
        theta[i]<-top -x[i]*beta
        # response is gamma
        y[i]~dgamma(shape, shape / (theta[i]))
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        shape ~ dlnorm(0,0.001) #dunif(0,1000)        
        top ~  dlnorm(0,0.001)

        # pearson residuals
        for (i in 1:N) {
         ExpY[i] <- theta[i] 
         VarY[i] <- shape/((shape/(shape / (theta[i])))^2)
         E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
         ysim[i] ~  dgamma(shape, shape / (theta[i]))
         Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
         D[i]    <- E[i]^2     #Squared residuals for original data
         Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = rgamma(1,0.2,0.001),
      shape = runif(1,0,10),
      top =rlnorm(1,log(quantile(mod.dat$y,probs = 0.25)),0.1))} 
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
        theta[i]<-top -x[i]*beta
        # response is gamma
        y[i]~dgamma(shape, shape / (theta[i]))
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        shape ~ dlnorm(0,0.001) #dunif(0,1000)
        top ~  dlnorm(0,0.001)

        # pearson residuals
        for (i in 1:N) {
         ExpY[i] <- theta[i] 
         VarY[i] <- shape/((shape/(shape / (theta[i])))^2)
         E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
         ysim[i] ~  dgamma(shape, shape / (theta[i]))
         Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
         D[i]    <- E[i]^2     #Squared residuals for original data
         Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = rgamma(1,0.2,0.001),
      shape = dlnorm(1,1/mean(mod.dat$y),1),
      top =rlnorm(1,log(quantile(mod.dat$y,probs = 0.25)),0.1))} 
  }
  
  # gamma y; gamma x -----
  if(x=="gamma" & y=="gamma"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top -x[i]*beta
        # response is gamma
        y[i]~dgamma(shape, shape / (theta[i]))
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        shape ~ dlnorm(0,0.001) #dunif(0,1000)
        top ~  dlnorm(0,0.001)

        # pearson residuals
        for (i in 1:N) {
         ExpY[i] <- theta[i] 
         VarY[i] <- shape/((shape/(shape / (theta[i])))^2)
         E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
         ysim[i] ~  dgamma(shape, shape / (theta[i]))
         Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
         D[i]    <- E[i]^2     #Squared residuals for original data
         Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])

        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = rgamma(1,0.2,0.001),
      shape = runif(1,0,10),
      top =rlnorm(1,log(quantile(mod.dat$y,probs = 0.25)),0.1))} 
  }
  
  # gaussian y; gamma x ----
 if(x=="gamma" & y=="gaussian"){
  sink("NECmod.txt")
  cat("
      model
      {
      
      # likelihood
      for (i in 1:N)
      {
      theta[i]<-top -x[i]*beta
      # response is gaussian
      
      y[i]~dnorm(theta[i],tau)
      }
      
      # specify model priors
      beta ~ dgamma(0.0001,0.0001)
      sigma ~ dunif(0, 20)  #sigma is the SD
      tau  = 1 / (sigma * sigma)  #tau is the reciprical of the variance 
      top ~  dnorm(0,0.1) # dnorm(0,0.001) #T(0,) 

      # pearson residuals
      for (i in 1:N) {
       ExpY[i] <- theta[i] 
       VarY[i] <- tau^2
       E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
      }

      # overdispersion
      for (i in 1:N) {
       ysim[i] ~  dnorm(theta[i],tau)
       Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
       D[i]    <- E[i]^2     #Squared residuals for original data
       Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
       }
       SS    <- sum(D[1:N])
       SSsim <- sum(Dsim[1:N])
      }
      ", fill=TRUE)
  sink()  #Make model in working directory
  
  init.fun <- function(mod.data=mod.data){list(
    beta = rgamma(1,0.2,0.001),
    sigma = runif(1, 0, 5),
    top =rnorm(1,min(mod.dat$y),1))}

 }

  # gaussian y; beta x ----
  if(x=="beta" & y=="gaussian"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top -x[i]*beta
        # response is gaussian
        
        y[i]~dnorm(theta[i],tau)
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        sigma ~ dunif(0, 20)  #sigma is the SD
        tau  = 1 / (sigma * sigma)  #tau is the reciprical of the variance 
        top ~  dnorm(0,0.1) # dnorm(0,0.001) #T(0,) 

      # pearson residuals
        for (i in 1:N) {
        ExpY[i] <- theta[i] 
        VarY[i] <- tau^2
        E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
        ysim[i] ~  dnorm(theta[i],tau)
        Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
        D[i]    <- E[i]^2     #Squared residuals for original data
        Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = rgamma(1,0.2,0.001),
      sigma = runif(1, 0, 5),
      top =rnorm(1,min(mod.dat$y),1))}
    
  }
  
  # gaussian y; gaussian x ----
  if(x=="gaussian" & y=="gaussian"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-top -x[i]*beta
        # response is gaussian
        
        y[i]~dnorm(theta[i],tau)
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        sigma ~ dunif(0, 20)  #sigma is the SD
        tau  = 1 / (sigma * sigma)  #tau is the reciprical of the variance 
        top ~  dnorm(0,0.1) # dnorm(0,0.001) #T(0,) 

      # pearson residuals
        for (i in 1:N) {
        ExpY[i] <- theta[i] 
        VarY[i] <- tau^2
        E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
        ysim[i] ~  dnorm(theta[i],tau)
        Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
        D[i]    <- E[i]^2     #Squared residuals for original data
        Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = rgamma(1,0.2,0.001),
      sigma = runif(1, 0, 5),
      top =rnorm(1,min(mod.dat$y),1))}
    
  }
  
  # beta y; beta x ----
  if(x=="beta" & y=="beta"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {

        # response is beta
        y[i]~dbeta(shape1[i], shape2[i])
        shape1[i] <- theta[i] * phi
        shape2[i]  <- (1-theta[i]) * phi
        theta[i]<-top -x[i]*beta
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        t0 ~ dnorm(0, 0.010)
        phi <- exp(t0)

        # pearson residuals
        for (i in 1:N) {
         ExpY[i] <- theta[i] 
         VarY[i] <- (shape1[i]*shape2[i])/((shape1[i]+shape2[i])^2*(shape1[i]+shape2[i]+1))
         E[i] <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
         ysim[i] ~ dbeta(shape1[i], shape2[i])
         Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
         D[i]    <- E[i]^2     #Squared residuals for original data
         Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = rgamma(1,0.2,0.001),
      t0 = rnorm(0,100),
      top =quantile(mod.dat$y, probs=0.25))} 
  }
  
  # beta y; gamma x ----
  if(x=="gamma" & y=="beta"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        
        # response is beta
        y[i] ~ dbeta(shape1[i], shape2[i])
        shape1[i] <- theta[i] * phi
        shape2[i]  <- (1-theta[i]) * phi
        theta[i]<-top -x[i]*beta
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        t0 ~ dnorm(0, 0.010)
        phi <- exp(t0)
        top ~  dunif(0.0001,0.99)

        # pearson residuals
        for (i in 1:N) {
         ExpY[i] <- theta[i] 
         VarY[i] <- (shape1[i]*shape2[i])/((shape1[i]+shape2[i])^2*(shape1[i]+shape2[i]+1))
         E[i] <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
         ysim[i] ~ dbeta(shape1[i], shape2[i])
         Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
         D[i]    <- E[i]^2     #Squared residuals for original data
         Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = rgamma(1,0.2,0.001),
      t0 = rnorm(0,100),
      top = quantile(mod.dat$y, probs=0.25))}
  }

  # beta y; gaussian x ----
  if(x=="gaussian" & y=="beta"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        
        # response is beta
        y[i]~dbeta(shape1[i], shape2[i])
        shape1[i] <- theta[i] * phi
        shape2[i]  <- (1-theta[i]) * phi
        theta[i]<-top -x[i]*beta
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        t0 ~ dnorm(0, 0.010)
        phi <- exp(t0)
        top ~  dunif(0.0001,0.99)

        # pearson residuals
        for (i in 1:N) {
        ExpY[i] <- theta[i] 
        VarY[i] <- (shape1[i]*shape2[i])/((shape1[i]+shape2[i])^2*(shape1[i]+shape2[i]+1))
        E[i] <- (y[i] - ExpY[i]) / sqrt(VarY[i])
  }
        
        # overdispersion
        for (i in 1:N) {
        ysim[i] ~ dbeta(shape1[i], shape2[i])
        Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
        D[i]    <- E[i]^2     #Squared residuals for original data
        Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = rgamma(1,0.2,0.001),
      t0 = rnorm(0,100),
      top = quantile(mod.dat$y, probs=0.25))}
  }
  
  # negbin y; gaussian x ----
  if(x=="gaussian" & y=="negbin"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-size/(size+ top -x[i]*beta)
        # response is begative binomial
        y[i]~dnegbin(theta[i], size)
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        size ~ dunif(0,50)
        top ~ dgamma(1,0.001)
        
        # pearson residuals
        for (i in 1:N) {
        ExpY[i] <- theta[i] 
        VarY[i] <- theta[i] + theta[i] * theta[i] / size
        E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
        ysim[i] ~  dnegbin(theta[i], size)
        Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
        D[i]    <- E[i]^2     #Squared residuals for original data
        Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])
        
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = runif(1,0.0001,0.999), #rlnorm(1,0,1), #rgamma(1,0.2,0.001),
      size=runif(1, 0.1, 40),
      top =rpois(1,min(mod.dat$y)))} #
  } 

    # negbin y; gamma x ----
  if(x=="gamma" & y=="negbin"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-size/(size + top -x[i]*beta)
        # response is begative binomial
        y[i]~dnegbin(theta[i], size)
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        size ~ dunif(0,50)
        top ~ dgamma(1,0.001)
        
        # pearson residuals
        for (i in 1:N) {
        ExpY[i] <- theta[i] 
        VarY[i] <- theta[i] + theta[i] * theta[i] / size
        E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
        ysim[i] ~  dnegbin(theta[i], size)
        Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
        D[i]    <- E[i]^2     #Squared residuals for original data
        Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])
        
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = runif(1,0.0001,0.999), #rlnorm(1,0,1), #rgamma(1,0.2,0.001),
      size=runif(1, 0.1, 40),
      top =rpois(1,min(mod.dat$y)))} #
  }  
 
  # negbin y; beta x ----
  if(x=="beta" & y=="negbin"){
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-size/(size + top -x[i]*beta)
        # response is begative binomial
        y[i]~dnegbin(theta[i], size)
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        size ~ dunif(0,50)
        top ~ dgamma(1,0.001)
        
        # pearson residuals
        for (i in 1:N) {
        ExpY[i] <- theta[i] 
        VarY[i] <- theta[i] + theta[i] * theta[i] / size
        E[i]    <- (y[i] - ExpY[i]) / sqrt(VarY[i])
        }
        
        # overdispersion
        for (i in 1:N) {
        ysim[i] ~  dnegbin(theta[i], size)
        Esim[i] <- (ysim[i] - ExpY[i]) / sqrt(VarY[i])
        D[i]    <- E[i]^2     #Squared residuals for original data
        Dsim[i] <- Esim[i]^2  #Squared residuals for simulated data
        }
        SS    <- sum(D[1:N])
        SSsim <- sum(Dsim[1:N])
        
        }
        ", fill=TRUE)
    sink()  #Make model in working directory
    
    init.fun <- function(mod.data=mod.data){list(
      beta = runif(1,0.0001,0.999),
      size=runif(1, 0.1, 40),
      top =rpois(1,min(mod.dat$y)))} #
  }  
  
  # return the initial function
  if(exists("init.fun")){
      return(init.fun) 
  }
  else{
    sbot(paste("jagsEC50 does not currently support ", x, " distributed concentration data with ", y, 
                 " distributed response data. Please check this is the correct distribution to use, and if so
          feel free to contact the developers to request to add this distribution", sep=""))
  }
}