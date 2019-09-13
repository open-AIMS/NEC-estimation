#' check.mixing
#'
#' Generates a plot of MCMC chains and ACF function for a Jags model output.
#'
#' @param  J a jag model fit as returned by a call to jags fit.jagsNEC
#'
#' @export
#' @return A plot of MCMC chains and ACF diagrames for each element in params.


check.mixing <- function(J){
  x <- J$BUGSoutput$sims.array
  params <- J$parameters.to.save
  num.chains <- ncol(x[,,params[1]])
  within.cv <- lapply(params, FUN=function(i){
     x1          <- as.vector(x[,,i])
     chain.id     <- rep(1:num.chains, each = nrow(x[,,i]))
     unlist(lapply(1:num.chains,FUN=function(k){
        sd(x1[which(chain.id==k)])/mean(x1[which(chain.id==k)])     
      }))       
    }) 
  names(within.cv) <- params
  
  between.cv <- lapply(params, FUN=function(i){
     x1          <- as.vector(x[,,i])
     chain.id     <- rep(1:num.chains, each = nrow(x[,,i]))
     sd(x1)/mean(x1) 
    }) 
  names(between.cv) <- params
  
  cv.test <- lapply(params, FUN=function(i){
    between.cv[[i]]>quantile(within.cv[[i]], probs = 0.75)
  })
  names(cv.test) <- params
    
  cv.ratio <- mean(unlist(between.cv)/unlist(lapply(within.cv, FUN=mean))) 
    
  return(list('within.cv'=within.cv,
              'between.cv'=between.cv,
              'cv.test'=cv.test,
              'cv.ratio'=cv.ratio))  
}

