
library(R2jags)
require(tidyverse)
source("R/check_chains.R")
source("R/check_mixing.R")
source("R/Write_jags_model.R")
source("R/Predict_fitted_vals.R")
source("R/Fit_jagsNEC.R")
source("R/plot_jagsNEC.R")
source("R/extract_ECx.R")
path <- "C:/Users/rfisher/OneDrive - Australian Institute of Marine Science/Documents/AIMS/EcologicalRiskModelling/Ecotoxicology/Ecotox_stats/CR-examples"


### Example from Gerards original NEC script ----
### #(https://github.com/gerard-ricardo/NECs/blob/master/NECs)
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

### test Heidi's necrosis example as a binomial (100 trials) -----

dat<-read.csv(paste(path,'results_CarterioAdultHCexp_noFWoutliers.csv',sep="/"), 
                          header=T, sep=',', strip.white=T) %>%
  mutate(Tr =Tr, #/100 deal with 0-1 x data later
         T96_p.healthy=100-T96_p.necrosed,
         Recov_p.healthy=100-Recov_p.necrosed,
         concentration=as.numeric(Tr),
         y.1=as.integer(floor(T96_p.healthy)),
         y.2=as.integer(floor(Recov_p.healthy)), #        
         trials=100)
View(dat)
str(dat)

out1 <- fit.jagsNEC(data=dat, 
                   x.var="concentration", 
                   y.var="y.1",
                   trials.var = "trials")


check.chains(out1) 
par(mfrow=c(1,1))
plot_jagsNEC(out1) 

out2 <- fit.jagsNEC(data=dat, 
                   x.var="concentration", 
                   y.var="y.2",
                   trials.var = "trials")

check.chains(out2)
par(mfrow=c(1,1))
plot_jagsNEC(out2) 


### test Heidi's necrosis example as a beta ----
dat<-read.csv(paste(path,'results_CarterioAdultHCexp_noFWoutliers.csv',sep="/"), 
              header=T, sep=',', strip.white=T) %>%
  mutate(Tr =Tr, # model x as beta
         T96_p.healthy=100-T96_p.necrosed,
         Recov_p.healthy=100-Recov_p.necrosed,
         concentration=as.numeric(Tr),
         y.1=T96_p.healthy/100,
         y.2=Recov_p.healthy/100)

# first for T96
out1 <- fit.jagsNEC(data=dat, 
                    x.var="concentration", 
                    y.var="y.1",
                    n.tries=2)


check.chains(out1) 
par(mfrow=c(1,1))
plot_jagsNEC(out1)

# second for Recovery
out2 <- fit.jagsNEC(data=dat, 
                    x.var="concentration", 
                    y.var="y.2",
                    n.tries=2)

check.chains(out2)
par(mfrow=c(1,1))
plot_jagsNEC(out2) 


### test Heidi's necrosis example as a binomial (trials as a function of area) ----
dat<-read.csv(paste(path,'results_CarterioAdultHCexp_noFWoutliers.csv',sep="/"), 
              header=T, sep=',', strip.white=T) %>%
  mutate(Tr=Tr/100,
         trials=round(TotalArea*10),
         T96_p.healthy=100-T96_p.necrosed,
         Recov_p.healthy=100-Recov_p.necrosed,
         concentration=as.numeric(Tr),
         y.1=as.integer(round(T96_p.healthy/100*trials)),
         y.2=as.integer(round(Recov_p.healthy/100*trials))
        )

# first for T96
out1 <- fit.jagsNEC(data=dat, 
                    x.var="concentration", 
                    y.var="y.1",
                    n.tries=2,
                    trials.var="trials")


check.chains(out1) 
par(mfrow=c(1,1))
plot_jagsNEC(out1)

# second for Recovery
out2 <- fit.jagsNEC(data=dat, 
                    x.var="concentration", 
                    y.var="y.2",
                    n.tries=2,
                    trials.var="trials")

check.chains(out2)
par(mfrow=c(1,1))
plot_jagsNEC(out2) 


#### test Beta using Gerard's example ----
prop.data <- read.table("https://pastebin.com/raw/123jq46d", header= TRUE,dec=",") %>%
  mutate(raw.x=log(as.numeric(as.character(raw.x))+1),
         resp=as.numeric(as.character(resp)))
out <- fit.jagsNEC(data=prop.data, 
                   x.var="raw.x", 
                   y.var="resp",
                   n.tries=1)

check.chains(out) 
par(mfrow=c(1,1))
plot_jagsNEC(out)

### test Heidi's necrosis example as a binomial for any necrosis at all (0,1) ----
dat<-read.csv(paste(path,'results_CarterioAdultHCexp_noFWoutliers.csv',sep="/"), 
              header=T, sep=',', strip.white=T) %>%
  mutate(Tr=Tr/100, # model x as beta
         T96_p.healthy=100-T96_p.necrosed,
         Recov_p.healthy=100-Recov_p.necrosed,
         concentration=as.numeric(Tr),
         y.1=as.integer(if_else(T96_p.healthy/100<1,0,1)),
         y.2=as.integer(if_else(Recov_p.healthy/100<1,0,1)),
         trials=1)

# first for T96
out1 <- fit.jagsNEC(data=dat, 
                    x.var="concentration", 
                    y.var="y.1",
                    trials.var="trials")

check.chains(out1) 
par(mfrow=c(1,1))
plot_jagsNEC(out1)

# second for Recovery
out2 <- fit.jagsNEC(data=dat, 
                    x.var="concentration", 
                    y.var="y.2",
                    trials.var="trials")

check.chains(out2)
par(mfrow=c(1,1))
plot_jagsNEC(out2, jitter.x=T) 

### Paul's sea urchins ----
binom.data <-  read.table(paste(path,"Data source R NEC and ECs sea urchin fertilization (Fisher, Ricardo, Fox).txt", sep="/"), header= TRUE,dec=",")
str(binom.data)
binom.data$raw.x <- as.numeric(as.character(binom.data$raw.x))
range(binom.data$raw.x)
par(mfrow=c(2,1))
hist(binom.data$raw.x)
hist(binom.data$suc/binom.data$tot)
out <- fit.jagsNEC(data=binom.data,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot")
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out, x.lab = "WAF (%)", y.lab = "Fertilization success (%)",log.x = "x")
extract_ECx(out)


### Gerard - 28/10/2019 example list ----
#hav4
data1 = read.table("https://pastebin.com/raw/v7Uy4iXs", header= TRUE,dec=",") %>% 
  mutate(raw.x=as.numeric(as.character(raw.x)))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot")
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out)
extract_ECx(out)

data1 = read.table("https://pastebin.com/raw/zfrUha88", header= TRUE,dec=",") %>%
mutate(raw.x=as.numeric(as.character(raw.x)))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot")
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out)
extract_ECx(out)


#paul
# note this example was a bit problematic (failed to have good chain mixing). Modelling as a % of WAF improved the outcome substantially.
data1 = read.table("https://pastebin.com/raw/dKVi6L3t", header= TRUE,dec=",") %>% 
mutate(raw.x=as.numeric(as.character(raw.x))/100)
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", n.tries=1)
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out, x.lab = "Proportion WAF")
extract_ECx(out)

#bent4  * #
data1 = read.table("https://pastebin.com/raw/dUMSAvYi", header= TRUE,dec=",") %>%
mutate(raw.x=as.numeric(as.character(raw.x)))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", n.tries=1)
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out)
extract_ECx(out)

#orp4 ## No NEC
# note this example was a bit problematic (failed to have good chain mixing). Modelling as a % of WAF improved the outcome substantially.
data1 = read.table("https://pastebin.com/raw/BaCAP3Sr", header= TRUE,dec=",") %>% 
mutate(raw.x=as.numeric(as.character(raw.x))/100)
out.1 <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", n.tries=1)
check.chains(out.1)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out.1)
extract_ECx(out.1)
# backtransform NEC
(out.1$NEC)*100

# run the same example on a log scale - works a bit better.
data1 = read.table("https://pastebin.com/raw/BaCAP3Sr", header= TRUE,dec=",") %>%
  mutate(raw.x=log(as.numeric(as.character(raw.x))+1))
out.2 <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", n.tries=1)
check.chains(out.2)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out.2)
extract_ECx(out.2)
# backtransform NEC
exp(out.2$NEC)-1

#flo1  * 
### would not converge - try on log scale
data1 = read.table("https://pastebin.com/raw/CgYVQp6f", header= TRUE,dec=",") %>% 
mutate(raw.x=log(as.numeric(as.character(raw.x))))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", n.tries=1)
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out)
extract_ECx(out)

#flo2 
## would not converge - try on log scale
data1 = read.table("https://pastebin.com/raw/2jsiuqP2", header= TRUE,dec=",") %>% 
mutate(raw.x=log(as.numeric(as.character(raw.x))))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", n.tries=1)
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out)
extract_ECx(out)

#flo3
### would not converge - try on log scale
data1 = read.table("https://pastebin.com/raw/my6xYhfG", header= TRUE,dec=",") %>% 
mutate(raw.x=log(as.numeric(as.character(raw.x))))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", n.tries=1)
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out)
extract_ECx(out)

#flo4
### would not converge - try on log scale
data1 = read.table("https://pastebin.com/raw/atPPnT8v", header= TRUE,dec=",") %>% 
mutate(raw.x=log(as.numeric(as.character(raw.x))))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", n.tries=1)
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out)
extract_ECx(out)

#flo5
### would not run at all - try on log scale. Added error catching to jagsNEC for when a model won't fit and there is a large concentration range in x
data1 = read.table("https://pastebin.com/raw/C7U2JcSS", header= TRUE,dec=",") %>% 
mutate(raw.x=log(as.numeric(as.character(raw.x))))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", n.tries=1)
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out)
extract_ECx(out)

#flo6
data1 = read.table("https://pastebin.com/raw/Fp0uLBNX", header= TRUE,dec=",") %>% 
mutate(raw.x=as.numeric(as.character(raw.x)))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", n.tries=1)
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out)
extract_ECx(out)



