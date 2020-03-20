

data=dat #data1#count.data#binom.data##
x.var="x"#"light.stress" # "log.x" #"concentration" #
y.var="y"#"range01.col"# "col.intensity" # "scaled.col"#"suc" #"count" #"response" #"suc" #
trials.var=NA#"tot"#
model="ECxLinear" #"NECHormesis"#"ECxWeibull1"#"ECxWeibull1"#"NEC3param"#"NECsigmoidal"#"ECx4param"  #

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
init.value.warning=TRUE

CI=TRUE
posterior.median=TRUE
median.model=FALSE
add.NEC=TRUE

precision=10000
ECx.val=10
type="absolute"
prob.vals=c(0.5, 0.025, 0.975)
xform=NA
posterior = FALSE


X=out.ma

legend.loc="topright"
add.EC10=FALSE

lxform=NA

jitter.x=FALSE
jitter.y=FALSE
ylab="response"
xlab="concentration"
xlim=NA
xticks=NA

