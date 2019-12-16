
#' extract_ECx
#'
#' Extracts the predicted ECx value as desired from a jagsNEC model fit obeject
#'
#' @param  X a jag model fit as returned by a call to jags from fit.jagsNEC
#' 
#' @param ECx.val the desired percentage effect value. This must be a value between 1 and 99 (for type = "relative" and "absolute"), defaults to 10.
#' 
#' @param type a character vector indicating if relative or absolute values for the ECx should be calculated. Takes values 
#' of "relative",  "absolute" (the default) or "direct". For the default model type fit by jagsNEC, "relative" is 
#' calculated as the percentage decrease from the maximum value of the response (top) to the minimum predicted value 
#' of the response, "absolute" is calculated as the the percentage decrease from the maximum value of the response (top) 
#' to 0 (or bot for a 4 parameter model fit),  and "direct" provides a direct estimate of the x value for a given y.
#' 
#' @param precision The number of unique x values over which to find ECx - large values will make the ECx estimate more precise.
#' 
#' @param posterior A logical value indicating if the full posterior sample of calculated ECx values should be returned instead of just the median and 95 credible intervals.
#' 
#' @param xform A function to apply to the returned estimated concentration values
#' 
#' @export
#' @return A vector containing the estimated ECx value, including upper and lower 95 percent Credible Interval bounds
#' 
#' @details Please note that the estimated ECx value is based on the equivalent percentage decrease from the range of the highest to the lowest estimate value across the range of the observed concentration (x) values. If the concentration response relationship is such that the full range of observed responses is not captured (ie a complete decline response at the highest level of exposure), the estimated ECx values may be lower than if the full concentration-response curve were available. Note this is therefore a conservative value.

extract_ECx <- function(X, ECx.val=10, precision=10000, posterior = FALSE, type="absolute", xform=NA){
  if(type!="direct"){
   if(ECx.val<1 | ECx.val>99){
    stop("Supplied ECx.val is not in the required range. Please supply a percentage value between 1 and 99.")
   }   
  }
 
  if(X$y.type=="gaussian" & max(grep("4param", X$model))!= 1  & type=="absolute"){
    stop("Absolute ECx values are not valid for a gaussian response variable unless a 4 parameter model is fit") 
  }

  label <- paste("EC",ECx.val,sep="_")
  
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

  if(type=="absolute" & X$model == "4param"){
     ECx.out <- apply(pred.vals$posterior, MARGIN=2, FUN=function(y){
     range.y <- range(y)
     ECx.y <- max(range.y)-diff(range.y)*(ECx.val/100)
     ECx.x <- x.vec[which.min(abs(y-ECx.y))]
     return(ECx.x)  
     })     
  }
  
  if(type=="absolute" & X$model != "4param"){
    ECx.out <- apply(pred.vals$posterior, MARGIN=2, FUN=function(y){
      range.y <- c(0, max(y))
      ECx.y <- max(range.y)-diff(range.y)*(ECx.val/100)
      ECx.x <- x.vec[which.min(abs(y-ECx.y))]
      return(ECx.x)  
    })     
  }
  
  if(type=="direct"){
    ECx.out <- apply(pred.vals$posterior, MARGIN=2, FUN=function(y){
      range.y <- range(y)
      ECx.y <- ECx.val
      ECx.x <- x.vec[which.min(abs(y-ECx.y))]
      return(ECx.x)
      
    }) 
  }
  
  ECx.estimate <- quantile(ECx.out, probs=c(0.5, 0.025, 0.975))
  names(ECx.estimate) <- c(label, paste(label, "lw", sep="_"), paste(label, "up", sep="_"))

  if(class(xform)=="function"){
    ECx.estimate <- xform(ECx.estimate)
    ECx.out <- xform(ECx.out)
  }

  
  if(posterior==FALSE){
    return(ECx.estimate)
  }else{
    return(ECx.out)}

  
}

