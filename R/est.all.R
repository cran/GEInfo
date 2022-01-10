est.all<-function(Eols,Wols,Yols,active.set,alpha0,a0,b0,kappa1,kappa2,xi,epsilon,S_G,S_GE)
{
  q <- ncol(Eols)-1 # Eols includes intercept and X
  p <- ncol(Wols) / (q+1)
  s<-0
  repeat{
    est.ac<-update.para.svd (Eols,Wols,Yols,active.set,alpha0,a0,b0,kappa1,kappa2,xi,epsilon,S_G,S_GE)
    aold<-a0
    alphaold<-alpha0
    bold<-b0
    a0<-est.ac$a[-1]
    alpha0<-est.ac$a[1]
    b0<-est.ac$b
    s<-s+1
    if(s>=5) break
  }
  active.set<-sort(union(active.set,which(b0!=0)))
  S_H <- unique(S_GE[,1])
  S_G_all <- sort(union(S_G,S_H))
  beta.ind <- (q+1)*(S_G_all-1)+1
  jind <- S_GE[,1]
  kind <- S_GE[,2]
  gamma.ind <- (q+1)*(jind-1) + 1 + kind
  info.set <- sort(c(beta.ind,gamma.ind))
  acs<-union(active.set, info.set)
  AS_c<-setdiff(seq(length(b0)), acs)
  est.asc<-update.para.svd.c (Eols,Wols,Yols,AS_c,alpha0,a0,b0,kappa1,kappa2,xi,epsilon,S_G,S_GE)
  bnew=est.asc
  active.set <- which(bnew!=0)
  return(list(a=a0,b=bnew,alpha=alpha0,active.set=active.set))
}
