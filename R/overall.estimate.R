overall.estimate  <- function (E,W,Y,active.set,alpha0,a0,b0,family,kappa1,kappa2,xi,epsilon,S_G,S_GE)
  {
  taylor <- taylor.exp(E,W,Y,alpha0,a0,b0,family)
  weight.IWLS <- taylor$weight
  sqrt.wt <-  as.vector (sqrt(weight.IWLS))
  Y.IWLS <- taylor$Ystar
  Yols <- sqrt.wt * Y.IWLS
  Eols <- sqrt.wt * cbind(1,E)
  Wols <- sqrt.wt * W
  est.new <- est.all(Eols,Wols,Yols,active.set,alpha0,a0,b0,kappa1,kappa2,xi,epsilon,S_G,S_GE)
  a0 <- est.new$a
  b0 <- est.new$b
  alpha0 <- est.new$alpha
  active.set <- est.new$active.set
  return(list(alphanew=alpha0,anew=a0,bnew=b0,active.set=active.set))
}
