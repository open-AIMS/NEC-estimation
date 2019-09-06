

data= count.data#binom.data
x.var="raw.x"
y.var= "count" #"suc"
#trials.var="tot"
x.type = "gamma"
y.type = "poisson" # #binomial"#"gaussian" #
params = c("top", "beta", "NEC")
burnin = 1000
n.iter = 2000



CI=TRUE
posterior.median=TRUE
median.model=FALSE
add.NEC=TRUE

X=out
