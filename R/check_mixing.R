#' check.mixing
#'
#' Calculates the within and between chain mixing and tests of the between is greater than within based on the supplied quantile probability value.
#'
#' @param  J a jag model fit as returned by a call to jags fit.jagsNEC
#'
#' @param  prob.val the probability value to pass to quantile to assess the between chain coefficient of variation against the within chain coefficient of variation. Please inspect the model code for more details.
#'
#' @export
#' @return A list of within and between coefficient of variations for parameter, the outcome of the test if between is greater than within, and the ratio ratio.


check.mixing <- function(J, prob.val=0.75){
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
    between.cv[[i]]>quantile(within.cv[[i]], probs = prob.val)
  })
  names(cv.test) <- params
    
  cv.ratio <- mean(unlist(between.cv)/unlist(lapply(within.cv, FUN=mean))) 
    
  return(list('within.cv'=within.cv,
              'between.cv'=between.cv,
              'cv.test'=cv.test,
              'cv.ratio'=cv.ratio))  
}

