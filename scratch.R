

data=dat #count.data#binom.data##
x.var="light.stress" #log.x" # "concentration" #
y.var="col.intensity" # "suc" #"count" #"response" #"suc" # 
#model="NECHormesis"#"NEC3param"#"NECsigmoidal"#"ECx4param"  #
trials.var=NA#"tot"
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
xform=NA

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
