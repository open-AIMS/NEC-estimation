
devtools::install_github("AIMS/NEC-estimation")
require(jagsNEC)

library(R2jags)
require(tidyverse)
require(readxl)
path <- "C:/Users/rfisher/OneDrive - Australian Institute of Marine Science/Documents/AIMS/EcologicalRiskModelling/Ecotoxicology/Ecotox_stats/CR-examples"


source("R/check_chains.R")
source("R/check_mixing.R")
source("R/Write_jags_NEC3paramMod.R")
source("R/Write_jags_NECsigmoidalMod.R")
source("R/Write_jags_NEC4paramMod.R")
source("R/Write_jags_ECx4paramMod.R")
source("R/Write_jags_NECHormesisMod.R")
source("R/Write_jags_ECxWeibull1Mod.R")
source("R/Write_jags_ECxWeibull2Mod.R")
source("R/Predict_fitted_vals.R")
source("R/Fit_jagsNEC.R")
source("R/Fit_jagsMANEC.R")
source("R/plot_jagsNEC.R")
source("R/plot_jagsNECfit.R")
source("R/extract_ECx.R")
source("R/wi.R")


#----test_dat1.csv ----
dat<-read.csv(paste(path,'test_dat1.csv',sep="/"))
out <- fit.jagsNEC(data=dat,
                   x.var="light.stress", 
                   y.var="col.intensity",
                   model="NECHormesis")#"ECxWeibull2") #

check.chains(out)

par(mfrow=c(1,1))
plot(out)

out.ma <- fit.jagsMANEC(data=dat, 
                          x.var="light.stress", 
                          y.var="col.intensity")

check.chains(out.ma)



#----test_dat2.csv-----

dat<-read.csv(paste(path,'test_dat2.csv',sep="/"))
out <- fit.jagsNEC(data=dat,
                   x.var="light.stress", 
                   y.var="scaled.col", n.tries=3,
                   #model="NEC3param")
                   model="NECHormesis")
check.chains(out)

par(mfrow=c(1,1))
plot(out)


out.ma <- fit.jagsMANEC(data=dat, 
                        x.var="light.stress", 
                        y.var="col.intensity")

#----test_dat3.csv-----

dat<-read.csv(paste(path,'test_dat3.csv',sep="/"))
out <- fit.jagsNEC(data=dat,
                   x.var="light.stress", 
                   y.var="range01.col", 
                   n.tries=1,
                   model="NECHormesis")
                 # model="NEC3param")

check.chains(out)
par(mfrow=c(1,1))
plot(out)


out.ma <- fit.jagsMANEC(data=dat,
                   x.var="light.stress",  n.tries=1,
                   y.var="range01.col")


par(mfrow=c(1,1))
plot(out.ma)

#----test_dat4.csv-----
dat<-read.csv(paste(path,'test_dat4.csv',sep="/"))
out3 <- fit.jagsNEC(data=dat,
                   x.var="light.stress", 
                   y.var="scaled.col", n.tries=1,
                   #model="NEC3param")
                   model="NECHormesis")
check.chains(out3)

par(mfrow=c(1,1))
plot(out3)

#-----test_dat5.csv-----
dat<-read.csv(paste(path,'test_dat5.csv',sep="/"))
out <- fit.jagsNEC(data=dat,
                   x.var="light.stress", 
                   y.var="scaled.col", n.tries=1,
                   #model="NEC3param")
                   model="NECHormesis")
check.chains(out)

par(mfrow=c(1,1))
plot(out)

#----test_dat6.csv-----
dat<-read.csv(paste(path,'test_dat6.csv',sep="/"))
out <- fit.jagsNEC(data=dat,
                   x.var="light.stress", 
                   y.var="scaled.col", n.tries=1,
                   #model="NEC3param")
                   model="NECHormesis")
check.chains(out)

par(mfrow=c(1,1))
plot(out)

### testing what to do with hormesis ----
data1 <- read.table("https://pastebin.com/raw/rF2biHYG", header= TRUE,dec=",") %>%
  mutate(raw.x = as.numeric(as.character(raw.x)),
         log.x = log(raw.x),
         suc = as.integer(as.character(suc)),
         tot = as.integer(as.character(tot)),
         prop = suc/tot) %>%
  drop_na() 

data1$obs <- factor(formatC(1:nrow(data1), flag="0", width = 3))# first need to make an observation row to soon randomise
plot(data1$raw.x, data1$suc/data1$tot, log = 'x')   #pooled by spp and specta (tenuis not updates

out <- fit.jagsNEC(data=data1,
                   x.var="log.x",  
                   y.var="suc",
                   model="ECxWeibull1",
                   #over.disp = TRUE,
                   trials.var = "tot")
check.chains(out)

par(mfrow=c(1,1))
plot(out)

out.ma <- fit.jagsMANEC(data=data1, 
                        x.var="log.x", 
                        y.var="suc",
                        trials.var = "tot",
                        over.disp = TRUE)

### Testing/troubleshooting alternative models ----
#source("R/Write_jags_NECHormesisMod.R")

data1 <-  read.table("https://pastebin.com/raw/dUMSAvYi", header= TRUE,dec=",")  %>%
  mutate(raw.x=as.numeric(as.character(raw.x)),
         log.x=log(raw.x))

out.ma <- fit.jagsMANEC(data=data1,
                    x.var="log.x",
                    y.var="suc",
                    over.disp = TRUE,
                    trials.var = "tot")

out1 <- fit.jagsNEC(data=data1,
                   x.var="log.x",
                   y.var="suc",
                   model="ECx4param",
                   over.disp = TRUE,
                   trials.var = "tot")

out2 <- fit.jagsNEC(data=data1,
                   x.var="log.x",
                   y.var="suc",
                   model="NEC4param",
                   over.disp = TRUE,
                   trials.var = "tot")
out3 <- fit.jagsNEC(data=data1,
                         x.var="log.x",
                         y.var="suc",
                         model="NECsigmoidal",
                         over.disp = TRUE,
                         trials.var = "tot", 
                         n.tries=20)
out4 <- fit.jagsNEC(data=data1,
                    x.var="log.x",
                    y.var="suc",
                    model="NECHormesis",
                    over.disp = TRUE,
                    trials.var = "tot")
out5 <- fit.jagsNEC(data=data1,
                    x.var="log.x",
                    y.var="suc",
                    model="NEC3param",
                    over.disp = TRUE,
                    trials.var = "tot")

pdf("compare_NECvals.pdf", height=8)
par(mfrow=c(3,2), mar=c(2,2,1,1), oma=c(2,2,0,0))
plot(out1)
abline(v=out1$NEC, lty=c(2,1,2), col="blue")
legend("bottomleft", legend="ECx4param", bty="n")
legend("right", legend="NEC (indirect)", col="blue", lty=1, bty="n")
plot(out5)
legend("bottomleft", legend="NEC3param", bty="n")
plot(out2)
legend("bottomleft", legend="NEC4param", bty="n")
plot(out3)
legend("bottomleft", legend="NECsigmoidal", bty="n")
plot(out4)
legend("bottomleft", legend="NECHormesis", bty="n")
mtext(text="Log (Concentration)", side=1, outer=TRUE)
mtext(text="Response", side=2, outer=TRUE)

dev.off()

out.list <- list(out1=out1, out2=out2, out3=out3,out4=out4,out5=out5)

require(custom.functions.pkg)

data.frame(DICw=unlist(round(wi(unlist(lapply(out.list, FUN=function(x){x$DIC}))),3)),
           pD=unlist(lapply(out.list, FUN=function(x){x$pD})))


out1$NEC
out2$NEC

check.chains(out)
par(mfrow=c(1,1))
plot(out)#, lxform = exp, xlim=c(-1,3), xticks = c(-1,0,1,2,3))
extract_ECx(out, type="direct", ECx.val=0.8)#, xform = exp)
extract_ECx(out, ECx.val=50)

### Testing/troubleshooting the Hockey Model v2

data1 <-  read.table("https://pastebin.com/raw/dKVi6L3t", header= TRUE,dec=",")  %>%
  mutate(raw.x=as.numeric(as.character(raw.x)),
         log.x=log(raw.x))
out.ma=fit.jagsMANEC(data=data1,
                     x.var="log.x",#"raw.x", # 
                     y.var="suc",
                     over.disp=1,
                     trials.var = "tot")
out <- fit.jagsNEC(data=data1,
                   x.var="log.x",#"raw.x", # 
                   y.var="suc",
                   model="NECHormesis",
                   trials.var = "tot")
check.chains(out)

par(mfrow=c(1,1))
plot(out)
plot(data1$log.x,out$residuals)


out1 <- fit.jagsNEC(data=data1,
                    x.var="log.x",
                    y.var="suc",
                    model="ECx4param",
                    over.disp = TRUE,
                    trials.var = "tot")

out2 <- fit.jagsNEC(data=data1,
                    x.var="log.x",
                    y.var="suc",
                    model="NEC4param",
                    over.disp = TRUE,
                    trials.var = "tot")
par(mfrow=c(1,2))
plot(out1)
plot(out2)
out1$NEC.p 
out2$NEC.p
out2$NEC


### Testing/troubleshooting alternative models
dat <- read.table(paste(path,'marie_test_hockey.txt',sep="/"), sep="\t", header=T) %>% 
  dplyr::select(raw.x, count) %>%
  mutate(SGR=(log(count)-log10(3353))/3,
         log.x=log(raw.x)) %>%
  na.omit()

out <- fit.jagsNEC(data=dat,
                   x.var="log.x",
                   y.var="SGR", 
                   y.type="gaussian",
                   model="NECHormesis",                   
                   n.tries=20)

check.chains(out)

par(mfrow=c(1,1))
plot(out, legend.loc = "bottomleft")
abline(v=out$NEC, col="orange", lty=c(2,1,2))

hist(out$sims.list$NEC)


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
                        trials.var="tot",
                        model="ECx4param" )

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
                   model="ECx4param",
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
                   model="ECx4param")
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
                   model="NECHormesis",
                   n.tries=20)

check.chains(out)
mean(out$sims.list$SS > out$sims.list$SSsim)

par(mfrow=c(1,1))
plot(out)


#### other testing ####################################
# now try on the logit scale just to experiment
binom.data$prop <- binom.data$suc/binom.data$tot
require(car)

binom.data$logit.prop <- logit(binom.data$prop)
out <- fit.jagsNEC(data=binom.data, 
                   x.var="log.x", 
                   y.var="logit.prop",
                   model="NECHormesis")
check.chains(out)

par(mfrow=c(1,1))
plot(out)

extract_ECx(out, ECx.val = 50)

### now test all Heidi's examples
all.files <- list.files(path)
files <- all.files[grep(".csv",all.files)]

pdf("testingECx4param.pdf",onefile = T)

for(f in 1:length(files)){
  dat <- read.csv(paste(path,files[f], sep="/"))
  out <- try(fit.jagsNEC(data=dat, 
                     x.var="concentration", 
                     y.var="response", 
                     model="ECx4param",
                     burnin=1000,
                     n.iter.update = 10000))
  #check.chains(out)
  if(class(out)!="try-error"){
      par(mfrow=c(1,1))
      plot(out)
      mtext(side=3,text=f,outer=T)
  }

  
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
                    model="ECx4param")

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
                    model="ECx4param")


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
                   model="ECx4param",
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
                   model="ECx4param")
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
                   model="NEC4param")
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
  mutate(raw.x=as.numeric(as.character(raw.x)),
         log.x=log(raw.x))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot")
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out)
extract_ECx(out)
out$over.disp

out <- fit.jagsNEC(data=data1,
                   x.var="log.x",
                   y.var="suc",
                   trials.var="tot",
                   model="NECHormesis")
out$over.disp
plot(out)

# try using beta
out <- fit.jagsNEC(data=data1,
                   x.var="log.x",
                   y.var="suc",
                   trials.var="tot",
                   over.disp = TRUE,
                   model="ECx4param")
out$over.disp
plot(out)



data1 = read.table("https://pastebin.com/raw/zfrUha88", header= TRUE,dec=",") %>%
mutate(raw.x=as.numeric(as.character(raw.x)),
       log.x=log(raw.x))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot")
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out)
out$over.disp
# overdispersed
out <- fit.jagsNEC(data=data1,
                   x.var="log.x",
                   y.var="suc",
                   trials.var="tot",
                   over.disp=T,
                   model="NEC4param") 
plot(out)
out$over.disp
extract_ECx(out)



#paul
# note this example was a bit problematic (failed to have good chain mixing). Modelling as a % of WAF improved the outcome substantially.
data1 = read.table("https://pastebin.com/raw/dKVi6L3t", header= TRUE,dec=",") %>% 
mutate(raw.x=as.numeric(as.character(raw.x))/100)
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", 
                   n.tries=1,
                   model="NEC4param")
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out, x.lab = "Proportion WAF")
out$over.disp
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
plot(out)
extract_ECx(out)
mean(out$sims.list$SS > out$sims.list$SSsim)

#orp4 ## No NEC
# note this example was a bit problematic (failed to have good chain mixing). Modelling as a % of WAF improved the outcome substantially.
data1 = read.table("https://pastebin.com/raw/BaCAP3Sr", header= TRUE,dec=",") %>% 
mutate(raw.x=as.numeric(as.character(raw.x))/100,
       log.x=log(raw.x*100+1))
out.1 <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", n.tries=1)
check.chains(out.1)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out.1)
extract_ECx(out.1)
# backtransform NEC
(out.1$NEC)*100
mean(out.1$sims.list$SS > out.1$sims.list$SSsim)

# run the same example on a log scale - works a bit better.
out.2 <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", 
                   n.tries=1, 
                   x.type="beta", 
                   model="NEC4param",
                   over.disp=1)
check.chains(out.2)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out.2)
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
plot(out)
extract_ECx(out)
out$over.disp

#flo2 
## would not converge - try on log scale
data1 = read.table("https://pastebin.com/raw/2jsiuqP2", header= TRUE,dec=",") %>% 
mutate(raw.x=log(as.numeric(as.character(raw.x))))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", 
                   n.tries=1,
                   over.disp=TRUE,
                   model="NEC4param")
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out)
extract_ECx(out)
out$over.disp

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
plot(out)
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
plot(out)
extract_ECx(out)
mean(out$sims.list$SS > out$sims.list$SSsim)

#flo5
### would not run at all - try on log scale. Added error catching to jagsNEC for when a model won't fit and there is a large concentration range in x
data1 = read.table("https://pastebin.com/raw/C7U2JcSS", header= TRUE,dec=",") %>% 
mutate(raw.x=log(as.numeric(as.character(raw.x))))
out <- fit.jagsNEC(data=data1,
                   x.var="raw.x",
                   y.var="suc",
                   trials.var="tot", 
                   n.tries=1)
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out)
extract_ECx(out)
out$over.disp

#flo6
data1 = read.table("https://pastebin.com/raw/Fp0uLBNX", header= TRUE,dec=",") %>% 
mutate(raw.x=as.numeric(as.character(raw.x)),
       log.x=log(raw.x))
out <- fit.jagsNEC(data=data1,
                   x.var="log.x",
                   y.var="suc",
                   trials.var="tot", 
                   n.tries=1,
                   over.disp = TRUE)
check.chains(out)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out)
extract_ECx(out)
out$over.disp

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
                   trials.var="tot")


check.chains(out1)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out1)
extract_ECx(out1, precision = 10000)
out1$over.disp

out2 <- fit.jagsNEC(data=dat, 
                    x.var="log.x", 
                    y.var="suc",
                    trials.var="tot")


check.chains(out2)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(out2)
extract_ECx(out2)
out2$over.disp


### Flo take 2 ----

dat.diuron <- read_excel(paste(path,"/zoox data_use for R.xlsx", sep=""), sheet="diuron") %>%
  data.frame() %>%
  na.omit %>%
  mutate(sqrt.x=sqrt(raw.x),
         log.x=log(raw.x))
  

dat.metribuzin <- read_excel(paste(path,"/zoox data_use for R.xlsx", sep=""), sheet="metribuzin") %>%
  data.frame() %>%
  na.omit() %>%
  mutate(sqrt.x=sqrt(raw.x),
         log.x=log(raw.x))

head(dat.diuron)
head(dat.metribuzin)

par(mfrow=c(3,2), mar=c(3,2,0,0), oma=c(2,2,1,1))
plot(dat.diuron$raw.x, dat.diuron$resp)
plot(dat.metribuzin$raw.x, dat.metribuzin$resp,  yaxt="n")

plot(dat.diuron$sqrt.x, dat.diuron$resp)
plot(dat.metribuzin$sqrt.x, dat.metribuzin$resp,  yaxt="n")

plot(dat.diuron$log.x, dat.diuron$resp)
plot(dat.metribuzin$log.x, dat.metribuzin$resp, yaxt="n")

out.diuron <- fit.jagsNEC(data=dat.diuron, 
                          x.var="log.x", 
                          y.var="resp")

out.metribuzin <- fit.jagsNEC(data=dat.metribuzin, 
                          x.var="log.x", 
                          y.var="resp")



dat.propazine <- read_excel(paste(path,"/zoox data_use for R.xlsx", sep=""), sheet = "propazine") %>%
  data.frame() %>%
  na.omit %>%
  mutate(sqrt.x=sqrt(raw.x),
         log.x=log(raw.x),
         scaled.x=as.vector(scale(raw.x)))

out <- fit.jagsMANEC(data=dat.propazine, 
                             x.var="sqrt.x", 
                             y.var="resp", 
                             y.type="gaussian")


par(mfrow=c(1,1))
plot(out)
out$mod.stats


out <- fit.jagsNEC(data=dat.propazine, 
                   x.var="sqrt.x", 
                   y.var="resp", 
                   y.type="gaussian",
                   model="NECHormesis")




