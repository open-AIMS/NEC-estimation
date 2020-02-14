

data=data1#dat #count.data#binom.data##
x.var="log.x" #"light.stress" # "concentration" #
y.var="suc" #"col.intensity" # "count" #"response" #"suc" #
model="ECxWeibull1"#"NECHormesis"#"ECxWeibull1"#"NEC3param"#"NECsigmoidal"#"ECx4param"  #
trials.var="tot"#NA#
x.type = NA #"gamma"
y.type = NA #"gaussian"#"poisson" # #binomial"#"gaussian" #
params=c("top", "beta", "NEC", "SS", "SSsim")#, "slope")
burnin = 1000
n.iter = 2000
n.iter.update = 5000
n.tries=1
over.disp=FALSE
model.set=c("NEC3param", "NEC4param", "NECsigmoidal", "NECHormesis", "ECx4param")
name=""

CI=TRUE
posterior.median=TRUE
median.model=FALSE
add.NEC=TRUE

precision=10000

type="absolute"
prob.vals=c(0.5, 0.025, 0.975)

posterior = FALSE


X=out

legend.loc="topright"
add.EC10=FALSE
xform=NA
lxform=NA

jitter.x=FALSE
jitter.y=FALSE
ylab="response"
xlab="concentration"
xlim=NA
xticks=NA
