
#' extract_ECx
#'
#' Extracts the predicted ECx value as desired from a jagsNEC model fit obeject
#'
#' @param  X a jag model fit as returned by a call to jags from fit.jagsNEC
#' 
#' @param ECx.val the desired percentage effect value. This must be a value between 1 and 99, defaults to 10.
#' 
#' @export
#' @return A vector containing the estimated ECx value, including upper and lower 95% Credible Interval bounds
#' @details Please note that the estimated ECx value is based on the equivalent percentage decrease from the range of the highest to the lowest estimate value across the range of the observed concentration (x) values. If the concentration response relationship is such that the full range of observed responses is not captured (ie a complete decline response at the highest level of exposure), the estimated ECx values may be lower than if the full concentration-response curve were available. Note this is therefore a conservative value.

extract_ECx <- function(X, ECx.val=10, precision=1000, posterior = FALSE){
  if(ECx.val<1 | ECx.val>99){
    stop("Supplied ECx.val is not in the required range. Please supply a percentage value between 1 and 99.")}  

  pred.vals <- predict_NECbugsmod(X, precision=precision)
  
  x.vec <- pred.vals$'x' 
  
  ECx.out <- apply(pred.vals$posterior, MARGIN=2, FUN=function(y){
    range.y <- range(y)
    ECx.y <- max(range.y)-diff(range.y)*(ECx.val/100)
    ECx.x <- x.vec[which.min(abs(y-ECx.y))]
    return(ECx.x)
  
  })
  
  label <- paste("EC",ECx.val,sep="")
  
  ECx.estimate <- quantile(ECx.out, probs=c(0.5, 0.025, 0.975))
  names(ECx.estimate) <- c(label, paste(label, "lw", sep="_"), paste(label, "up", sep="_"))
  if(posterior==FALSE){
    return(ECx.estimate)
  }else{  return(list(ECx.estimate=ECx.estimate, ECx.posterior=ECx.out))}

  
}


