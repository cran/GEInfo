overall.est.final <- function (E,W,Y,active.set,alpha0,a0,b0,family,kappa1,kappa2,xi,epsilon,S_G,S_GE) {
  taylor <- taylor.exp(E,W,Y,alpha0,a0,b0,family)
  weight.IWLS <- taylor$weight
  sqrt.wt <-  as.vector (sqrt(weight.IWLS))
  Y.IWLS <- taylor$Ystar
   Yols <- sqrt.wt * Y.IWLS
  Eols <- sqrt.wt * cbind(1,E)
  Wols <- sqrt.wt * W
  s<-0
  repeat{
    est.ac<-update.para.svd (Eols,Wols,Yols,active.set,alpha0,a0,b0,kappa1,kappa2,xi,epsilon,S_G,S_GE)
    aold<-round(a0,4)
    alphaold<-round(alpha0,4)
    bold<-round(b0,4)
    a0<-round(est.ac$a[-1],4)
    alpha0<-round(est.ac$a[1],4)
    b0<-round(est.ac$b,4)
    theta.new <- c(a0,b0)
    theta.old <- c(aold,bold)
    value <- round(abs(theta.old-theta.new)/(1+abs(theta.old)),4)
    s=s+1
    if (max(value) <= 1e-3 | s>=100 ) break
  }
  return(list(alphanew=alpha0,anew=a0,bnew=b0))
}
