#' plot.jagsNEC
#'
#' Generates a plot of a fitted jags NEC model, as returned by fit.jagsNEC.
#' 
#' @param X the jags NEC model fit (as returned by fit.jagsNEC)
#' 
#' @param CI a logical value indicating if confidence intervals on the model fit should be plotted, calculated as the upper and lower bounds of the individual predicted values from all posterior samples
#'
#' @param posterior.median a logical value indicating if the posterior median of the model fit should be plotted, calculated as the median of the individual predicted values from all posterior samples
#'
#' @param median.model a logical value indicating of the fitted model caluclated from the median estimates of the NEC, top and beta parameters should be plotted. This is the fitted model as shown in Fox 2010.
#'
#' @param add.NEC a logocal value indicating if the estimated NEC values and 95% credible intervals should be added to the plot.
#'
#' @export
#' @return a plot of the fitted model

plot.jagsNEC <- function(X,
  CI=TRUE,
  posterior.median=TRUE,
  median.model=FALSE,
  add.NEC=TRUE){
 
  y.type <- out$y.type
  
  if(y.type=="binomial"){
     plot(X$mod.dat$x,X$mod.dat$y/X$mod.dat$trials, ylab="response", pch=16, 
       col=adjustcolor(1, alpha=0.25), cex=1.5, xlab="concentration") 
  }else{
    plot(X$mod.dat$x,X$mod.dat$y, ylab="response", pch=16, 
       col=adjustcolor(1, alpha=0.25), cex=1.5, xlab="concentration")     
  }

  
  abline(v=X$NEC, col = "red", lty=c(3,2,3))   
  
  if(CI==TRUE){
    lines(X$pred.vals$x, X$pred.vals$up, lty=2) 
    lines(X$pred.vals$x, X$pred.vals$lw, lty=2)  
  }
  if(posterior.median==TRUE){
    lines(X$pred.vals$x, X$pred.vals$y)
  }
  if(median.model==TRUE){
    lines(X$pred.vals$x, X$pred.vals$y.m, col="red")  
  }
  if(add.NEC==TRUE){
        legend("topright", bty="n",
           legend=paste("NEC: ", signif(X$NEC[2],3), 
                        " (", signif(X$NEC[1],3),"-", signif(X$NEC[3],3),")",sep=""))
  }

  
}
