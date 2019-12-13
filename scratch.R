

data=binom.data#data1 #count.data#
x.var="raw.x" # "concentration" #
y.var="suc" #"response" #"count" #"suc" # 
model=""#basic4param"  #"Hockey"
trials.var="tot"
x.type = NA #"gamma"
y.type = "gaussian"#"poisson" # #binomial"#"gaussian" #
params=c("top", "beta", "NEC", "SS", "SSsim")
burnin = 1000
n.iter = 2000
n.iter.update = 5000
n.tries=10

name=""

CI=TRUE
posterior.median=TRUE
median.model=FALSE
add.NEC=TRUE



over.disp=FALSE

X=out

legend.loc="topright"
add.EC10=FALSE
xform=NA
lxform=NA

jitter.x=FALSE
jitter.y=FALSE
ylab="response"
xlab="concentration"