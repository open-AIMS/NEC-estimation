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
#' @param add.NEC a logical value indicating if the estimated NEC values and 95\% credible intervals should be added to the plot.
#'
#' @param xform a function to be applied as a transformation of the x data.
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

plot_jagsNEC <- function(X,  CI=TRUE,  posterior.median=TRUE,  median.model=FALSE,  
                         add.NEC=TRUE, xform=NA,
                         jitter.x=FALSE, jitter.y=FALSE, 
                         y.lab="response", 
                         x.lab="concentration", log.x="", ...){
  
  # check if y.type is binomial
  y.type <- X$y.type
  if(y.type=="binomial"){
    y.dat <- X$mod.dat$y/X$mod.dat$trials}else{
    y.dat <- X$mod.dat$y}
  
  # check if a transformation is required for x
  if(class(xform)=="function"){
    x.dat <- xform(X$mod.dat$x)
    NEC <- xform(X$NEC)
    x.vec <- xform(X$pred.vals$x)
  }else{
    x.dat <- X$mod.dat$x
    NEC <- X$NEC
    x.vec <- X$pred.vals$x
  }
  
  if(jitter.x==TRUE){x.dat <- jitter(x.dat)}
  if(jitter.y==TRUE){y.dat <- jitter(y.dat)}  
  
  plot(x.dat, y.dat, ylab=y.lab, pch=16, 
       col=adjustcolor(1, alpha=0.25), cex=1.5, xlab=x.lab, log=log.x)   
  
  abline(v=NEC, col = "red", lty=c(3,1,3))   
  
  if(CI==TRUE){
    lines(x.vec, X$pred.vals$up, lty=2) 
    lines(x.vec, X$pred.vals$lw, lty=2)  
  }
  if(posterior.median==TRUE){
    lines(x.vec, X$pred.vals$y)
  }
  if(median.model==TRUE){
    lines(x.vec, X$pred.vals$y.m, col="red")  
  }
  if(add.NEC==TRUE){
        legend("topright", bty="n",
           legend=paste("NEC: ", signif(NEC[2],2), 
                        " (", signif(NEC[1],2),"-", signif(NEC[3],2),")",sep=""))
  }
  
  
  warning("The function plot_jagsNEC has been replaced with a generic plot function. You can now just call plot()")

}