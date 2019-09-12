

data=dat #count.data#binom.data
x.var="concentration" #"raw.x"
y.var="response" #"count" #"suc"
trials.var=NA#"tot"
x.type = NA #"gamma"
y.type = NA #"poisson" # #binomial"#"gaussian" #
params = c("top", "beta", "NEC")
burnin = 1000
n.iter = 2000



CI=TRUE
posterior.median=TRUE
median.model=FALSE
add.NEC=TRUE

X=out

