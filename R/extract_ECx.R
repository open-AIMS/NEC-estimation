
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

extract_ECx <- function(X, ECx.val=10){
  if(ECx.val<1 | ECx.val>99){
    stop("Supplied ECx.val is not in the required range. Please supply a percentage value between 1 and 99.")}  

  x.vec <- X$pred.vals$'x' 
  y <- X$pred.vals$y
  up <- X$pred.vals$up
  lw <- X$pred.vals$lw
  
  range.y <- range(y)
  range.up <- range(up)
  range.lw <- range(lw)
  
  ECx.y <- max(range.y)-diff(range.y)*(ECx.val/100)
  ECx.x <- x.vec[which.min(abs(y-ECx.y))]
  
  ECx.y.up <- max(range.up)-diff(range.up)*(ECx.val/100)
  ECx.x.up <- x.vec[which.min(abs(up-ECx.y.up))]
 
  ECx.y.lw <- max(range.lw)-diff(range.lw)*(ECx.val/100)
  ECx.x.lw <- x.vec[which.min(abs(lw-ECx.y.lw))] 
  label <- paste("EC",ECx.val,sep="")
  
  ECx.estimate <- c(ECx.x,ECx.x.lw,ECx.x.up)
  names(ECx.estimate) <- c(label, paste(label, "lw", sep="_"), paste(label, "up", sep="_"))
  return(ECx.estimate)
  
}


