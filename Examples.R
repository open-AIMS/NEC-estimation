
library(R2jags)
source("R/check_chains.R")
source("R/Write_jags_model.R")
source("R/Predict_fitted_vals.R")
source("R/Fit_jagsNEC.R")
source("R/plot_jagsNEC.R")

### Example from Gerards original NEC script (https://github.com/gerard-ricardo/NECs/blob/master/NECs)
binom.data <-  read.table("https://pastebin.com/raw/zfrUha88", header= TRUE,dec=",")
binom.data$raw.x <- as.numeric(as.character(binom.data$raw.x))
range(binom.data$raw.x)
hist(binom.data$raw.x)
# for x, lowest concentration is 0.1, highest is 400, right skewed and on the continuous scale. Model x as gamma (current default).
# for y, this is clearly binomial, with totals and successes
out <- fit.jagsNEC(data=binom.data, 
                        x.var="raw.x", 
                        y.var="suc", 
                        trials.var="tot",  
                        y.type = "binomial",
                        burnin=5000)

check.chains(out)

par(mfrow=c(1,1))
plot.jagsNEC(out)



# now try the log and see if it makes any difference
binom.data$prop <- binom.data$suc/binom.data$tot











min.x <- min(mod.dat$x)
max.x <- max(mod.dat$x)
x.seq <- seq(min.x, max.x, length=100)
 



mod.dat <<- list(
  x = binom.data$raw.x,   # concentration
  y = binom.data$suc, # response (successes)
  trials = binom.data$tot, # number of "trials"
  N = nrow(binom.data))  # Sample size

params <- c("top", "beta", "NEC")


init.fun <- write.jags.NECmod(x="gamma",y="binomial")


J1 <- jags(data       = mod.dat,
           inits      = init.fun,
           parameters = params,
           model      = "NECmod.txt",
           n.thin     = 10,
           n.chains   = 5,
           n.burnin   = 5000,
           n.iter     = 6000)

pred.vals <- predict.NECbugsmod(out)


NEC <-  quantile(out$sims.list$NEC,c(0.025, 0.5, 0.975))
top <-  quantile(out$sims.list$top,c(0.025, 0.5, 0.975))
beta <-  quantile(out$sims.list$beta,c(0.025, 0.5, 0.975)) 
y.pred.m <- predict.NECmod(x.vec=x.seq, NEC=NEC["50%"], top=top["50%"], beta=beta["50%"]) 

plot(out$mod.dat$x,out$mod.dat$y/out$mod.dat$trials, ylab="response", pch=16, 
     col=adjustcolor(1, alpha=0.25), cex=1.5, xlab="concentration")  
abline(v=NEC, col = "red", lty=c(3,2,3))   

if(CI==TRUE){
 lines(out$pred.vals$x, out$pred.vals$up, lty=2) 
 lines(out$pred.vals$x, out$pred.vals$lw, lty=2)  
}
if(posterior.median==TRUE){
 lines(out$pred.vals$x, out$pred.vals$y)
}
if(median.model==TRUE){
 lines(out$pred.vals$x, out$pred.vals$y.m, col="red")  
}
if(add.NEC==TRUE){
 legend("topright", bty="n",
       legend=paste("NEC: ", round(NEC[2],2), 
                    " (", round(NEC[1],2),"-", round(NEC[3],2),")",sep=""))
}
