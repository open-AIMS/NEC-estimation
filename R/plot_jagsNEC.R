#' plot.jagsNEC
#'
#' Generates a plot of a fitted jags NEC model, as returned by fit.jagsNEC.
#' 
#' @param X the jags NEC model fit as returned by fit.jagsNEC.
#' 
#' @param CI a logical value indicating if confidence intervals on the model fit should be plotted, calculated as the upper and lower bounds of the individual predicted values from all posterior samples
#'
#' @param posterior.median a logical value indicating if the posterior median of the model fit should be plotted, calculated as the median of the individual predicted values from all posterior samples
#'
#' @param median.model a logical value indicating if the fitted model calculated from the median estimates of the NEC, top and beta parameters should be plotted. This is the fitted model as shown in Fox 2010.
#'
#' @param add.NEC a logocal value indicating if the estimated NEC values and 95\% credible intervals should be added to the plot.
#'
#' @param x.jitter a logical value indicating if the x data points on the plot should be jittered.
#'
#' @param y.jitter a logical value indicating if the y data points on the plot should be jittered.
#' 
#' @param y.lab A character string containing a label for the y-axis label
#' 
#' @param x.lab A character string containing a label for the x-axis
#' 
#' @param log.x A character string indicating if the x-axis should be plotted on a log scale. Leave as the default for a plot without a log x-axis, or "x" will provide a plot with the x-axis on the log scale.
#'
#' @export
#' @return a plot of the fitted model

plot_jagsNEC <- function(X,  CI=TRUE,  posterior.median=TRUE,  median.model=FALSE,  add.NEC=TRUE, 
                         jitter.x=FALSE, jitter.y=FALSE, y.lab="response", x.lab="concentration", log.x=""){
 
  y.type <- X$y.type
  
  if(jitter.x==TRUE){X$mod.dat$x <- jitter(X$mod.dat$x)}
  if(jitter.y==TRUE){X$mod.dat$y <- jitter(X$mod.dat$y)}
  
  if(y.type=="binomial"){
     plot(X$mod.dat$x,X$mod.dat$y/X$mod.dat$trials, ylab=y.lab, pch=16, 
       col=adjustcolor(1, alpha=0.25), cex=1.5, xlab=x.lab, log=log.x) 
  }else{
    plot(X$mod.dat$x,X$mod.dat$y, ylab=y.lab, pch=16, 
       col=adjustcolor(1, alpha=0.25), cex=1.5, xlab=x.lab, log=log.x)     
  }

  abline(v=X$NEC, col = "red", lty=c(3,1,3))   
  
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