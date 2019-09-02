




predict.NECmod <- function(x.vec, NEC, top, beta){

  x.seq.pre <-  x.vec[which(x.vec<=NEC)] #seq(min.x, NEC.m, length=20)
  x.seq.post <- x.vec[which(x.vec>NEC)] # seq(NEC.m, max.x, length=20)
  y.pred.pre <- (rep(top, length(x.seq.pre)))
  y.pred.post <- (top*exp(-beta*(x.seq.post-NEC)))
  y.pred <- c(y.pred.pre, y.pred.post)
  return(y.pred)

}

predict.NECbugsmod <- function(mod.dat, X){
  min.x <- min(mod.dat$x)
  max.x <- max(mod.dat$x)
  x.seq <- seq(min.x, max.x, length=100)

  pred.vals.out <- do.call("cbind",lapply(1:X$n.sims,FUN=function(x){
    predict.NECmod(x.vec=x.seq,
                   NEC=X$sims.list$NEC[x],
                   top=X$sims.list$top[x],
                   beta=X$sims.list$beta[x])}))
  m.vals <- apply(pred.vals.out, MARGIN=1, FUN=quantile, probs=0.5)
  up.vals <- apply(pred.vals.out, MARGIN=1, FUN=quantile, probs=0.975)
  lw.vals <- apply(pred.vals.out, MARGIN=1, FUN=quantile, probs=0.025)

  return(list(
    x=x.seq,
    y=m.vals,
    up=up.vals,
    lw=lw.vals
  ))
}
