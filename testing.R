
#devtools::install_github("AIMS/NEC-estimation")
#require(jagsNEC)
library(R2jags)
require(tidyverse)
source("R/check_chains.R")
source("R/check_mixing.R")
source("R/Write_jags_model.R")
source("R/Write_jags_Hockey_model.R")
source("R/Write_jags_4param_model.R")
source("R/Write_jags_basic4param_model_2.R")
source("R/Predict_fitted_vals.R")
source("R/Fit_jagsNEC.R")
source("R/plot_jagsNEC.R")
source("R/plot_jagsNECfit.R")
source("R/extract_ECx.R")
path <- "C:/Users/rfisher/OneDrive - Australian Institute of Marine Science/Documents/AIMS/EcologicalRiskModelling/Ecotoxicology/Ecotox_stats/CR-examples"

### Testing/troubleshooting alternative models

data1 <-  read.table("https://pastebin.com/raw/dUMSAvYi", header= TRUE,dec=",")  %>%
  mutate(raw.x=as.numeric(as.character(raw.x)),
         log.x=log(raw.x))

out <- fit.jagsNEC(data=data1,
                   x.var="log.x",
                   y.var="suc",
                   trials.var = "tot",
                   model="basic4param",                   
                   n.tries=2,
                   over.disp=T)
check.chains(out)

par(mfrow=c(1,1))
plot(out)#, lxform = exp)
extract_ECx(out, type="direct", ECx.val=0.8)#, xform = exp)
extract_ECx(out, xform = exp, ECx.val=50)

### Testing/troubleshooting the Hockey Model v2

data1 <-  read.table("https://pastebin.com/raw/dKVi6L3t", header= TRUE,dec=",")  %>%
  mutate(raw.x=as.numeric(as.character(raw.x)),
         log.x=log(raw.x))

out <- fit.jagsNEC(data=data1,
                   x.var="log.x",#"raw.x", 
                   y.var="suc",
                   trials.var = "tot",
                   model="basic4param",                   
                   n.tries=5)
check.chains(out)

par(mfrow=c(1,1))
plot(out)
plot(data1$log.x,out$residuals)

### Testing/troubleshooting the Hockey Model v3
dat <- read.table(paste(path,'marie_test_hockey.txt',sep="/"), sep="\t", header=T) %>% 
  dplyr::select(raw.x, count) %>%
  mutate(SGR=(log(count)-log10(3353))/3,
         log.x=log(raw.x)) %>%
  na.omit()

out <- fit.jagsNEC(data=dat,
                   x.var="log.x",
                   y.var="SGR", y.type="gaussian",
                   model="basic4param",                   
                   n.tries=20)

check.chains(out)

par(mfrow=c(1,1))
plot(out, legend.loc = "bottomleft")




### Example from Gerards original NEC script ----
### #(https://github.com/gerard-ricardo/NECs/blob/master/NECs)
#1) NEC - binomial data  (a count out of a total, think %survival of individuals, % settlement)
binom.data <-  read.table("https://pastebin.com/raw/zfrUha88", header= TRUE,dec=",")
binom.data$raw.x <- as.numeric(as.character(binom.data$raw.x))
binom.data$log.x <- log(binom.data$raw.x)
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
plot(out)

out$over.disp
mean(out$sims.list$SS > out$sims.list$SSsim)
extract_ECx(out, type="relative")

# try log x
out <- fit.jagsNEC(data=binom.data, 
                   x.var="log.x", 
                   y.var="suc", 
                   model="basic4param",
                   trials.var="tot", over.disp=T)

check.chains(out)

par(mfrow=c(1,1))
plot(out, lxform=exp)

out$over.disp

extract_ECx(out, type="absolute")

####################################################
#2) NEC - count data  (think invidivuals, cells etc
count.data = read.table("https://pastebin.com/raw/ENgNSgf7", header= TRUE,dec=",") %>%
  mutate(raw.x=as.numeric(as.character(raw.x)),
         log.x=log(raw.x))
str(count.data)

range(count.data$raw.x)
par(mfrow=c(2,1))
hist(count.data$raw.x)
hist(count.data$count)
out <- fit.jagsNEC(data=count.data, 
                   x.var="log.x", 
                   y.var="count", 
                   model="basic4param")
check.chains(out)
mean(out$sims.list$SS > out$sims.list$SSsim)

par(mfrow=c(1,1))
plot(out)
extract_ECx(out, type="relative")
extract_ECx(out)
            
###################################################
#3) NEC - measured/continuous data  (think anything on the metric scale)
measure.data = read.table("https://pastebin.com/raw/pWeS6x0n", header= TRUE,dec=",") %>%
 mutate(raw.x = as.numeric(as.character(raw.x)),
        measure = as.numeric(as.character(measure)),
        log.x=log(raw.x))

out <- fit.jagsNEC(data=measure.data, 
                   x.var="log.x", 
                   y.var="measure", 
                   model="basic4param",
                   n.tries=1, prob.val=0.05)

check.chains(out)
mean(out$sims.list$SS > out$sims.list$SSsim)

par(mfrow=c(1,1))
plot_jagsNEC(out)


#### other testing ####################################
# now try on the logit scale just to experiment
binom.data$prop <- binom.data$suc/binom.data$tot
require(car)

binom.data$logit.prop <- logit(binom.data$prop)
out <- fit.jagsNEC(data=binom.data, 
                   x.var="log.x", 
                   y.var="logit.prop",
                   model="basic4param")
check.chains(out)

par(mfrow=c(1,1))
plot(out)

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
check.chains(out)
mean(out1$sims.list$SS > out1$sims.list$SSsim)


check.chains(out1) 
par(mfrow=c(1,1))
plot(out1)

# second for Recovery
out2 <- fit.jagsNEC(data=dat, 
                    x.var="concentration", 
                    y.var="y.2",
                    n.tries=2,
                    model="basic4param")

check.chains(out2)
par(mfrow=c(1,1))
plot(out2) 
mean(out2$sims.list$SS > out2$sims.list$SSsim)

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
                    trials.var="trials",
                    model="basic4param")


check.chains(out1) 
par(mfrow=c(1,1))
plot(out1)

# second for Recovery
out2 <- fit.jagsNEC(data=dat, 
                    x.var="concentration", 
                    y.var="y.2",
                    n.tries=2,
                    trials.var="trials")

check.chains(out2)
par(mfrow=c(1,1))
plot_jagsNEC(out2) 


### test Beta using Gerard's example ----
prop.data <- read.table("https://pastebin.com/raw/123jq46d", header= TRUE,dec=",") %>%
  mutate(raw.x=log(as.numeric(as.character(raw.x))+1),
         resp=as.numeric(as.character(resp)))
out <- fit.jagsNEC(data=prop.data, 
                   x.var="raw.x", 
                   y.var="resp",
                   model="basic4param",
                   n.tries=1)

check.chains(out) 
par(mfrow=c(1,1))
plot(out)
mean(out$sims.list$SS > out$sims.list$SSsim)

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
mean(out1$sims.list$SS > out2$sims.list$SSsim)

# second for Recovery
out2 <- fit.jagsNEC(data=dat, 
                    x.var="concentration", 
                    y.var="y.2",
                    trials.var="trials", 
                    over.disp=TRUE,
                    model="Hockey",
                    n.tries = 20)

check.chains(out2)
par(mfrow=c(1,1))
plot(out2, jitter.x=T) 
mean(out2$sims.list$SS > out2$sims.list$SSsim)

### Paul's sea urchins ----
binom.data <-  read.table(paste(path,"Data source R NEC and ECs sea urchin fertilization (Fisher, Ricardo, Fox).txt", 
                                sep="/"), header= TRUE,dec=",") %>%
  mutate(raw.x= as.numeric(as.character(raw.x)),
         log.x=log(raw.x),
         prop = suc/tot)

str(binom.data)
binom.data$
range(binom.data$raw.x)
par(mfrow=c(2,1))
hist(binom.data$raw.x)
hist(binom.data$suc/binom.data$tot)
out <- fit.jagsNEC(data=binom.data,
                   x.var="log.x",
                   y.var="suc",
                   trials.var="tot")
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out, x.lab = "WAF (%)", y.lab = "Fertilization success (%)", lxform=exp)
extract_ECx(out)
out$over.disp
out$summary
par(mfrow=c(1,1), mar=c(4,4,1,1))


# Paul's sea urchins as beta
out <- fit.jagsNEC(data=binom.data,
                   x.var="log.x",
                   y.var="prop", 
                   n.tries=2,
                   model="basic4param")
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out, x.lab = "WAF (%)", y.lab = "Fertilization success (%)")
extract_ECx(out)
out$over.disp
out$summary
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(binom.data$raw.x, out$residuals)

# Paul's sea urchins as beta using over.disp=TRUE
out <- fit.jagsNEC(data=binom.data,
                   x.var="log.x",
                   y.var="suc",
                   trials.var="tot", 
                   over.disp=T,
                   model="4param")
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out, x.lab = "WAF (%)", y.lab = "Fertilization success (%)",log.x = "x")
extract_ECx(out)
out$over.disp
out$summary
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(binom.data$raw.x, out$residuals)


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
mean(out$sims.list$SS > out$sims.list$SSsim)

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
mean(out$sims.list$SS > out$sims.list$SSsim)

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
mean(out$sims.list$SS > out$sims.list$SSsim)

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
mean(out$sims.list$SS > out$sims.list$SSsim)

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
mean(out.1$sims.list$SS > out.1$sims.list$SSsim)

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
mean(out.2$sims.list$SS > out.2$sims.list$SSsim)

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
mean(out$sims.list$SS > out$sims.list$SSsim)

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
mean(out$sims.list$SS > out$sims.list$SSsim)

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
mean(out$sims.list$SS > out$sims.list$SSsim)

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
mean(out$sims.list$SS > out$sims.list$SSsim)

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
mean(out$sims.list$SS > out$sims.list$SSsim)

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
mean(out$sims.list$SS > out$sims.list$SSsim)

### Marie pesticide MS ------
t0.dat <- unlist(list(
  Diuron= 3440,
  Metribuzin= 3320,
  Hexazinone= 3353,
  Bromacil= 3120,
  Tebuthiuron= 3971,
  Simazine= 2837,
  Propazine= 3114,
  Imazapic= 4966,
  Haloxyfop= 7255,
  H2.4.D= 3067))

t0.dat

path <- "C:/Users/rfisher/OneDrive - Australian Institute of Marine Science/Documents/AIMS/EcologicalRiskModelling/Ecotoxicology/Marie_Thomas/"
all.files <- list.files(path)
dat.files <- all.files[grep("_Rs.txt", all.files)]
names.vec <- gsub("2,4-D", "H2.4.D", gsub("_growth_Rs.txt","",dat.files))
dat.list <- list()
pdf(file="Marie_herbicides.pdf", onefile = T)
for(d in 1:length(dat.files)){
  t0.d <- t0.dat[names.vec[d]]
  dat.d <- read.table(paste(path, dat.files[d], sep=""), sep="\t", header=T) %>% 
    dplyr::select(raw.x, count) %>%
    mutate(SGR=(log(count)-log(t0.d))/3,
           log.x=log(raw.x))
  
  out <- try(fit.jagsNEC(data=na.omit(dat.d),
                     x.var="log.x",
                     y.var="SGR"), silent=T)
  if(class(out)!="try-error"){
     check.chains(out)
  par(mfrow=c(1,1), mar=c(4,4,1,1))
  plot_jagsNEC(out, x.lab = names.vec[d], y.lab = "Specific growth rate") 
  }else{
    plot(dat.d$log.x, dat.d$SGR, xlab=names.vec[d], "Specific growth rate")
  }
 dat.list <- c(dat.list, list(dat.d))
}
dev.off()
names(dat.list) <- names.vec

### Flo Diazinon ----
dat <- read.table(paste(path,"NECEstimation_Diazinon.txt", sep="/"), header = T) %>%
  mutate(log.x=log(raw.x))
out1 <- fit.jagsNEC(data=dat, 
                   x.var="raw.x", 
                   y.var="suc",
                   trials.var="tot",
                   burnin=100000,
                   n.iter.update = 10000)


check.chains(out1)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out1)
extract_ECx(out1, precision = 10000)
mean(out1$sims.list$SS > out1$sims.list$SSsim)

out2 <- fit.jagsNEC(data=dat, 
                    x.var="log.x", 
                    y.var="suc",
                    trials.var="tot",
                    burnin=10000,
                    n.iter.update = 1000)


check.chains(out2)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot_jagsNEC(out2)
extract_ECx(out2)
mean(out2$sims.list$SS > out2$sims.list$SSsim)
