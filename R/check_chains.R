#' check.chains
#'
#' Generates a plot of MCMC chains and ACF function for a Jags model output.
#'
#' @param  X a jag model fit as returned by a call to jags from R2jags
#'
#' @param  params A vector of names indicating the parameters that were traced during the jags fit. For the NEC jags model this is typically 'NEC','top' and 'beta', whicha are the defaults. 
#'
#' @param name and optional character string indicating the label to be placed at the top of the plotting window
#'
#' @export
#' @return A plot of MCMC chains and ACF diagrames for each element in params.

check.chains <- function(X, params = c("top", "beta", "NEC"), name=""){
  x <- X$sims.array
  num.chains <- ncol(x[,,params[1]])

  par(mfrow=c(length(params),2),mar=c(0,5,0.5,0.5),oma=c(4,0,2,0))
  for (i in 1:length(params)){
   	x1          <- as.vector(x[,,params[i]])
  	chain.id     <- rep(1:num.chains, each = nrow(x[,,params[i]]))
    num.lags <- length(acf(x[,,params[i]][,1], plot = FALSE)$lag)

    # plot the chains
    plot(1:nrow(x[,,i]),rep(NA,nrow(x[,,params[i]])), xaxt="n",
         ylim=range(x1), main="", xlab="", ylab=params[i])
    if(i == length(params)){axis(side=1)}
    for(k in 1:num.chains){
        lines(1:nrow(x[,,params[i]]), x1[which(chain.id==k)],col=k)}

    # plot the acf
    plot(1:num.lags,rep(NA,num.lags), xaxt="n",
        ylim=c(0,1),xlab="lag",ylab="correlation",main="")
    if(i == length(params)){axis(side=1)}
     for (j in 1:num.chains){
      acf.j <- acf(x[,,params[i]][,j], plot = FALSE)
      lines(acf.j$lag, acf.j$acf,col=j )
     }
    }
 mtext(name,side=3,outer=T)
}

