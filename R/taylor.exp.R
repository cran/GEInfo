taylor.exp<-function(E,W,Y,alpha0,a0,b0,family){
  n<-length(Y)
  weight<-rep(1,n)
  if(family == "gaussian") {
    weight <- weight
    Ystar <- Y
  }
   if(family=='poisson') {
    eta <- E%*%a0+W%*%b0+alpha0
    mu <- exp(eta)
    weight<-mu
    Ystar <- eta +(Y-mu) / weight
  }
  if(family=='binomial') {
    eta<-E%*%a0+W%*%b0+alpha0
    prob<-1/(1+exp(-eta))
    weight<-1/4
    mod.Y<-(Y-prob)/(weight)
    Ystar<-E%*%a0+W%*%b0+alpha0+mod.Y
  }
  return(list(weight=weight,Ystar=Ystar))
}
