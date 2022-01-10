
CGEInfo1 <- function(E,G,Y,family,lam1,lam2,
                    xi=6,epsilon=0,max.it=100,thresh=1e-03,S_G=NULL,S_GE=NULL)
{
  nlambda1<-length(lam1)
  nlambda2<-length(lam2)
  W<-matW(E,G)
  q<-ncol(E)
  p<-ncol(G)
  array.a <- array(0,dim=c(nlambda1,nlambda2,q))
  array.b <- array(0,dim=c(nlambda1,nlambda2,(q+1)*p))
  array.beta <- array(0,dim=c(nlambda1,nlambda2,p))
  array.gamma<-array(0,dim=c(nlambda1,nlambda2,p,q))
  array.est<-array(0,dim=c(nlambda1,nlambda2,(q+1)*(p+1)))
  mat.alpha <- matrix(0,nrow=nlambda1,ncol=nlambda2)
  if(family=='poisson')
  {
    fit0<-glmnet::cv.glmnet(x=cbind(E,W),y=Y, family=family)
    a0<-round(stats::coef(fit0)[2:(q+1)],3)
    b0<-round(stats::coef(fit0)[-(1:(q+1))],3)
    alpha0<-round(stats::coef(fit0)[1],3)
  } else {
    a0 <- rep(0,q)
    b0 <-rep(0,p*(1+q))
    alpha0 <- 0
  }
  active.set<-which(b0!=0)
  for(i in nlambda1:1) {
    for(j in nlambda2:1) {
      s <- 0
      repeat{
        est.new <- overall.estimate (E,W,Y,active.set,alpha0,a0,b0,family,kappa1=lam1[i],kappa2=lam2[j],xi,epsilon,S_G,S_GE)
        aold<-a0
        bold=b0
        alphaold=alpha0
        active.setold=active.set
        a0 <- round(est.new$anew,4)
        b0 <-round(est.new$bnew,4)
        alpha0<- round(est.new$alphanew,4)
        active.set <- est.new$active.set
        s=s+1
        if(s>=100 | identical(active.set,active.setold)) break
      }
      est.final<-overall.est.final(E,W,Y,active.set,alpha0,a0,b0,family,kappa1=lam1[i],kappa2=lam2[j],xi,epsilon,S_G,S_GE)
      a0=est.final$anew
      b0=est.final$bnew
      alpha0=est.final$alphanew
      array.a[i,j,] <- a0
      array.b[i,j,] <- b0
      mat.alpha[i,j] <- alpha0
      array.est[i,j,]<-c(alpha0,a0,b0)
      array.beta[i,j,]<-b0[(1:((q+1)*p))%%(q+1)==1]
      array.gamma [i,j,,]<- matrix(b0[(1:((q+1)*p))%%(q+1)!=1],ncol=q,nrow=p,byrow = T)
    }
  }
  est.res<-(list(a=array.a,beta=array.beta, gamma=array.gamma,b=array.b,alpha=mat.alpha,est.all=array.est,lambda1=lam1,lambda2=lam2))
  return(est.res)

}
