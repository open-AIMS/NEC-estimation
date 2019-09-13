
library(R2jags)
source("R/check_chains.R")
source("R/check_mixing.R")
source("R/Write_jags_model.R")
source("R/Predict_fitted_vals.R")
source("R/Fit_jagsNEC.R")
source("R/plot_jagsNEC.R")

### Example from Gerards original NEC script (https://github.com/gerard-ricardo/NECs/blob/master/NECs)
#1) NEC - binomial data  (a count out of a total, think %survival of individuals, % settlement)
binom.data <-  read.table("https://pastebin.com/raw/zfrUha88", header= TRUE,dec=",")
binom.data$raw.x <- as.numeric(as.character(binom.data$raw.x))
range(binom.data$raw.x)
hist(binom.data$raw.x)
# for x, lowest concentration is 0.1, highest is 400, right skewed and on the continuous scale. Model x as gamma (current default).
# for y, this is clearly binomial, with totals and successes
out <- fit.jagsNEC(data=binom.data, 
                        x.var="raw.x", 
                        y.var="suc", 
                        trials.var="tot")

check.chains(out)

par(mfrow=c(1,1))
plot_jagsNEC(out)


####################################################
#2) NEC - count data  (think invidivuals, cells etc
count.data = read.table("https://pastebin.com/raw/ENgNSgf7", header= TRUE,dec=",")
str(count.data)

count.data$raw.x <- as.numeric(as.character(count.data$raw.x))

range(count.data$raw.x)
par(mfrow=c(2,1))
hist(count.data$raw.x)
hist(count.data$count)
out <- fit.jagsNEC(data=count.data, 
                   x.var="raw.x", 
                   y.var="count")
check.chains(out)

par(mfrow=c(1,1))
plot_jagsNEC(out)

###################################################
#3) NEC - measured/continuous data  (think anything on the metric scale)
measure.data = read.table("https://pastebin.com/raw/pWeS6x0n", header= TRUE,dec=",")
measure.data
measure.data$raw.x <- as.numeric(as.character(measure.data$raw.x))
measure.data$measure <- as.numeric(as.character(measure.data$measure))

out <- fit.jagsNEC(data=measure.data, 
                   x.var="raw.x", 
                   y.var="measure", n.tries=1, prob.val=0.05)
check.chains(out)

par(mfrow=c(1,1))
plot_jagsNEC(out)


#### other testing ####################################
# now try on the logit scale just to experiment
binom.data$prop <- binom.data$suc/binom.data$tot
require(car)

binom.data$logit.prop <- logit(binom.data$prop)
out <- fit.jagsNEC(data=binom.data, 
                   x.var="raw.x", 
                   y.var="logit.prop")
check.chains(out)

par(mfrow=c(1,1))
plot_jagsNEC(out)

### now test all Heidi's examples
path <- "C:/Users/rfisher/OneDrive - Australian Institute of Marine Science/Documents/AIMS/EcologicalRiskModelling/Ecotoxicology/Ecotox_stats/CR-examples"
all.files <- list.files(path)
files <- all.files[grep(".csv",all.files)]

pdf("testing.pdf",onefile = T)

for(f in 1:length(files)){
  dat <- read.csv(paste(path,files[f], sep="/"))
  out <- fit.jagsNEC(data=dat, 
                     x.var="concentration", 
                     y.var="response", 
                     burnin=1000,
                     n.iter.update = 10000)
  check.chains(out)
  
  par(mfrow=c(1,1))
  plot_jagsNEC(out)
  mtext(side=3,text=f,outer=T)
  
}
dev.off()



