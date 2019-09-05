
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
require(car)

params <- c("top", "beta", "NEC", "alpha", "sigma")
binom.data$logit.prop <- logit(binom.data$prop)
out <- fit.jagsNEC(data=binom.data, 
                   x.var="raw.x", 
                   y.var="logit.prop", 
                   #trials.var="tot",  
                   y.type = "gaussian",
                   burnin=10000,
                   params=params)

check.chains(out, params=params)

par(mfrow=c(1,1))
plot.jagsNEC(out)






