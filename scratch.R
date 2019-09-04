

data=binom.data
x.var="raw.x"
y.var="suc"
trials.var="tot"
x.type = "gamma"
y.type = "binomial"
params = c("top", "beta", "NEC")
burnin = 1000
n.iter = 2000



CI=TRUE
posterior.median=TRUE
median.model=FALSE
add.NEC=TRUE


