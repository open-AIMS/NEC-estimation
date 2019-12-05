
#' extract_ECx
#'
#' Extracts the predicted ECx value as desired from a jagsNEC model fit obeject
#'
#' @param  X a jag model fit as returned by a call to jags from fit.jagsNEC
#' 
#' @param ECx.val the desired percentage effect value. This must be a value between 1 and 99, defaults to 10.
#' 
#' @param type a character vector indicating if relative or absolute values for the ECx should be calculated. Takes either "relative" or "absolute"
#' 
#' @param precision The number of unique x values over which to find ECx - large values will make the ECx estimate more precise.
#' 
#' @param posterior A logical value indicating if the full posterior sample of calculated ECx values should be returned instead of just the median and 95 credible intervals.
#' 
#' @export
#' @return A vector containing the estimated ECx value, including upper and lower 95 percent Credible Interval bounds
#' 
#' @details Please note that the estimated ECx value is based on the equivalent percentage decrease from the range of the highest to the lowest estimate value across the range of the observed concentration (x) values. If the concentration response relationship is such that the full range of observed responses is not captured (ie a complete decline response at the highest level of exposure), the estimated ECx values may be lower than if the full concentration-response curve were available. Note this is therefore a conservative value.

extract_ECx <- function(X, ECx.val=10, precision=10000, posterior = FALSE, type="relative"){
  if(ECx.val<1 | ECx.val>99){
    stop("Supplied ECx.val is not in the required range. Please supply a percentage value between 1 and 99.")
  } 
  if(is.na(match(X$y.type, c("binomial", "beta"))) & type=="absolute"){
    stop("Absolute ECx values are only valid for response variables that are contained from 0 to 1, 
         such as a binomial or beta") 
  }

  label <- paste("EC",ECx.val,sep="")
  
  pred.vals <- predict_NECbugsmod(X, precision=precision)
  
  x.vec <- pred.vals$'x' 
  
  if(type=="relative"){
   ECx.out <- apply(pred.vals$posterior, MARGIN=2, FUN=function(y){
     range.y <- range(y)
     ECx.y <- max(range.y)-diff(range.y)*(ECx.val/100)
     ECx.x <- x.vec[which.min(abs(y-ECx.y))]
     return(ECx.x)
  
   })    
  }

  if(type=="absolute"){
    ECx.out <- apply(pred.vals$posterior, MARGIN=2, FUN=function(y){
      range.y <- range(y)
      ECx.y <- 1-(ECx.val/100)
      ECx.x <- x.vec[which.min(abs(y-ECx.y))]
      return(ECx.x)
      
    })    
  }
  
  ECx.estimate <- quantile(ECx.out, probs=c(0.5, 0.025, 0.975))
  names(ECx.estimate) <- c(label, paste(label, "lw", sep="_"), paste(label, "up", sep="_"))
  if(posterior==FALSE){
    return(ECx.estimate)
  }else{
    return(ECx.out)}

  
}

