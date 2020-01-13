

data=data1 #count.data#binom.data##
x.var="log.x" #log.x" # "concentration" #
y.var="resp" # "suc" #"count" #"response" #"suc" # 
model="NEC3param"#"NECHormesis"#"NECsigmoidal"#"ECx4param"  #
trials.var=NA#"tot"
x.type = NA #"gamma"
y.type = NA #"gaussian"#"poisson" # #binomial"#"gaussian" #
params=c("top", "beta", "NEC", "SS", "SSsim", "slope")
burnin = 1000
n.iter = 2000
n.iter.update = 5000
n.tries=1
over.disp=FALSE

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
