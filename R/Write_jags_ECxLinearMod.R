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

write.jags.ECxLinear.mod <- function(x = "gamma", y, mod.dat) {

  # binomial y ----
  if (y == "binomial") {
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        logit(theta[i]) <- top - x[i]*beta
        #theta[i]<- top * exp(-beta*x[i])
        # response is binomial
        y[i]~dbin(theta[i],trials[i])
        }
        
        # specify model priors
        beta ~ dlnorm(0,0.1)
        top ~  dnorm(0,0.01)
         
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
        ", fill = TRUE)
    sink() # Make model in working directory

    init.fun <- function(mod.data = mod.data) {
      list(
        beta = rlnorm(1, 0, 1), #
        top = rnorm(1, 0, 1)
      )
    } # rgamma(1,0.2,0.001))}
    attributes(init.fun) <- list(link = "logit")
    return(init.fun)
  }

  # poisson y ----
  if (y == "poisson") {
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        log(theta[i])<- top - x[i]*beta
        #theta[i]<- top * exp(-beta*x[i])
        # response is poisson
        y[i] ~ dpois(theta[i])
        }
        
        # specify model priors
        beta ~ dlnorm(0,0.1)
        top ~  dnorm(0,0.001)
        
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
        ", fill = TRUE)
    sink() # Make model in working directory

    init.fun <- function(mod.data = mod.data) {
      list(
        beta = rlnorm(1, 0, 1),
        top = rnorm(1, 0, 1)
      )
    }
    attributes(init.fun) <- list(link = "log")
    return(init.fun)
  }

  # gamma y  -----
  if (y == "gamma") {
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        log(theta[i])<- top - x[i]*beta
        #theta[i]<- top * exp(-beta*x[i])
        # response is gamma
        y[i]~dgamma(shape, shape / (theta[i]))
        }
        
        # specify model priors
        beta ~ dlnorm(0,0.1)
        shape ~ dlnorm(0,0.1) #dunif(0,1000)        
        top ~  dnorm(0,0.001)

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
        ", fill = TRUE)
    sink() # Make model in working directory

    init.fun <- function(mod.data = mod.data) {
      list(
        beta = rlnorm(1, 0, 1),
        shape = runif(1, 0, 10),
        top = rnorm(1, 0, 1)
      )
    }
    attributes(init.fun) <- list(link = "log")
    return(init.fun)
  }

  # gaussian  ----
  if (y == "gaussian") {
    sink("NECmod.txt")
    cat("
      model
      {
      
      # likelihood
      for (i in 1:N)
      {
        theta[i]<- top - x[i]*beta
        #theta[i]<- top * exp(-beta*x[i])
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
      ", fill = TRUE)
    sink() # Make model in working directory

    init.fun <- function(mod.data = mod.data) {
      list(
        beta = rlnorm(1, 0, 1),
        sigma = runif(1, 0, 5),
        top = rnorm(1, 0, 1)
      )
    }
    attributes(init.fun) <- list(link = "identity")
    return(init.fun)
  }

  # beta y ------
  if (y == "beta") {
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
        logit(theta[i])<- top - x[i]*beta
        #theta[i]<- top * exp(-beta*x[i])
        }
        
        # specify model priors
        beta ~ dlnorm(0,0.1)
        t0 ~ dnorm(0, 0.010)
        phi <- exp(t0)
        top ~  dnorm(0,0.01)

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
        ", fill = TRUE)
    sink() # Make model in working directory

    init.fun <- function(mod.data = mod.data) {
      list(
        beta = rlnorm(1, 0, 1),
        t0 = rnorm(1, 0, 1),
        top = rnorm(1, 0, 1)
      )
    }
    attributes(init.fun) <- list(link = "logit")
    return(init.fun)
  }

  # negbin y ----
  if (y == "negbin") {
    sink("NECmod.txt")
    cat("
        model
        {
        
        # likelihood
        for (i in 1:N)
        {
        theta[i]<-size/(size + top - beta*x[i])

        # response is negative binomial
        y[i]~dnegbin(theta[i], size)
        }
        
        # specify model priors
        beta ~ dgamma(0.0001,0.0001)
        size ~ dunif(0,50)
        top ~ dgamma(1,0.001) # 
        
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
        ", fill = TRUE)
    sink() # Make model in working directory

    init.fun <- function(mod.data = mod.data) {
      list(
        beta = rgamma(1, 0.2, 0.001),
        size = runif(1, 0.1, 40),
        top = rlnorm(1, 0, 1)
      )
    }
    attributes(init.fun) <- list(link = "identity")
    return(init.fun)
  }
}
