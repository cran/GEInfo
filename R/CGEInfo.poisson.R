CGEInfo.poisson <- function(E,G,Y,family='poisson',lam1,lam2,
                            xi=6,epsilon=0,max.it=500,thresh=1e-03,S_G=NULL,S_GE=NULL)
{
  nlambda1<-length(lam1)
  nlambda2<-length(lam2)
  W<-matW(E,G)
  q<-ncol(E)
  p<-ncol(G)
  array.a <- array(0,dim=c(nlambda1,nlambda2,q))
  array.b <- array(0,dim=c(nlambda1,nlambda2,(q+1)*p))
  array.est<-array(0,dim=c(nlambda1,nlambda2,(q+1)*(p+1)))
  mat.alpha <- matrix(0,nrow=nlambda1,ncol=nlambda2)
  array.lam<-array(0,dim=c(nlambda1,nlambda2,2))
  fit0<-glmnet::cv.glmnet(x=cbind(E,W),y=Y, family=family)
  for(i in nlambda1:1)
  {
    for(j in nlambda2:1)
    {
      s <- 0
      a0<-round(stats::coef(fit0)[2:(q+1)],4)
      b0<-round(stats::coef(fit0)[-(1:(q+1))],4)
      alpha0<-round(stats::coef(fit0)[1],4)
      active.set<-which(b0!=0)
      repeat{
        try({est.new <-
          overall.estimate(E,W,Y,active.set,alpha0,a0,b0,family,kappa1=lam1[i],kappa2=lam2[j],xi,epsilon,S_G,S_GE)
        aold<-a0
        bold=b0
        alphaold=alpha0
        active.setold=active.set
        a0 <- round(est.new$anew,4)
        b0 <-round(est.new$bnew,4)
        alpha0<- round(est.new$alphanew,4)
        active.set <- est.new$active.set },silent=T)
        s=s+1
        if(s>=100 | identical(active.set,active.setold)) break
      }
       array.a[i,j,] <- a0
      array.b[i,j,] <- b0
      mat.alpha[i,j] <- alpha0
      array.lam[i,j,]<-c(lam1[i],lam2[j])
      array.est[i,j,]<-c(alpha0,a0,b0)
    }
  }
  est.res<-(list(a=array.a,b=array.b,alpha=mat.alpha,est.all=array.est,lambda1=lam1,lambda2=lam2))
  return(est.res)
}
