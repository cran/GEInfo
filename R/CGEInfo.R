#' @title  CGEInfo and GEsgMCP approaches with fixed tunings
#' @description Realize to estimate CGEInfo and GEsgMCP approaches at fixed tunings.
#' @param E Observed matrix of E variables, of dimensions n x q.
#' @param G Observed matrix of G variables, of dimensions n x p.
#' @param Y Response variable of length n. Quantitative for family="gaussian", or family="poisson" (non-negative count). For family="binomial" should be a factor with two levels.
#' @param family Model type: one of ("gaussian", "binomial", "poisson").
#' @param lam1 A user supplied lambda1.
#' @param lam2 A user supplied lambda2.
#' @param xi Tuning parameter of MCP penalty. Default is 6.
#' @param epsilon Tuning parameter of Ridge penalty which shrinks the coefficients of variables having prior information. Default is 0.
#' @param max.it Maximum number of iterations (total across entire path). Default is 500.
#' @param thresh Convergence threshold for group coordinate descent algorithm. The algorithm iterates until the change for each coefficient is less than thresh. Default is 1e-3.
#' @param S_G A user supplied vector, denoting the subscript of G variables which have prior information. Default is NULL.
#' @param S_GE A user supplied matrix, denoting the subscript of G-E interactions which have prior information.
#'             The first and second columns of S_GE represent the subscript of G variable and the subscript of E variable, respectively.
#'             For example, S_GE = matrix( c(1, 2), ncol = 2), which indicates that the 1st G and the 2nd E variables have an interaction effect on Y.
#'             Default is NULL. If both S_G and S_GE are NULL, no prior information is incorporated in the model, in which case function CGEInfo realizes GEsgMCP approach.
#'
#' @return An object of class "GEInfo" is returned, which is a list including the estimation results at fixed tunings.
#' \item{a}{Coefficient vector of length q for E variables.}
#' \item{b}{Coefficient vector of length (q+1)p for W (G variables and G-E interactions).}
#' \item{beta}{Coefficient vector of length p for G variables.}
#' \item{gamma}{Coefficient matrix of dimensions p*q for G-E interactions.}
#' \item{alpha}{Intercept.}
#' \item{coef}{A coefficient vector of length (q+1)*(p+1), including the estimates for \eqn{\alpha} (intercept), \eqn{a} (coefficients for all E variables), and \eqn{b} (coefficients for all G variables and G-E interactions).}
#' @references Wang X, Xu Y, and Ma S. (2019). Identifying gene-environment interactions incorporating prior information. Statistics in medicine, 38(9): 1620-1633. \doi{10.1002/sim.8064}

#' @export
#' @examples
#' n <- 30; p <- 5; q <- 2
#' E <- MASS::mvrnorm(n, rep(0,q), diag(q))
#' G <- MASS::mvrnorm(n, rep(0,p), diag(p))
#' W <- matW(E, G)
#' alpha <- 0; a <- seq(0.4, 0.6, length=q);
#' beta <- c(seq(0.2, 0.5, length=3),rep(0, p-3))  # coefficients of G variables
#' vector.gamma <- c(0.8, 0.5, 0, 0)
#' gamma <- matrix(c(vector.gamma, rep(0, p*q - length(vector.gamma))), nrow=p, byrow=TRUE)
#' mat.b.gamma <- cbind(beta, gamma)
#' b <- as.vector (t(mat.b.gamma))              # coefficients of G and G-E interactions
#' Y <- alpha + E %*% a + W %*% b + rnorm (n, 0, 0.5)
#' S_G <- c(1)
#' S_GE <- cbind(c(1), c(1))
#' fit1 <- CGEInfo(E, G, Y,family='gaussian', S_G=S_G, S_GE=S_GE,lam1=0.4,lam2=0.4)

CGEInfo <- function(E,G,Y,family,lam1,lam2,
                    xi=6,epsilon=0,max.it=500,thresh=1e-03,S_G=NULL,S_GE=NULL)
{
  W<-matW(E,G)
  q<-ncol(E)
  p<-ncol(G)
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
   s <- 0
      repeat{
        est.new <- overall.estimate (E,W,Y,active.set,alpha0,a0,b0,family,kappa1=lam1,kappa2=lam2,xi,epsilon,S_G,S_GE)
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
      est.final<-overall.est.final(E,W,Y,active.set,alpha0,a0,b0,family,kappa1=lam1,kappa2=lam2,xi,epsilon,S_G,S_GE)


      a0=est.final$anew
      if(is.null(colnames(E))) {Enames<-paste0('E',1:q)} else { Enames<-colnames(E)}
      if(is.null(colnames(G))) {Gnames<-paste0('G',1:p)} else { Gnames<-colnames(G)}
      names(a0)<-Enames

      b0=est.final$bnew
      index.beta0<- seq(length(b0))%%(q+1)==1
      names(b0)[index.beta0]<-Gnames
      names(b0)[!index.beta0]<-paste0(rep(Gnames,each=q),'-',rep(Enames,p))

      alpha0=est.final$alphanew
      beta0<-b0[index.beta0]
      gamma0<-matrix(b0[seq(length(b0))%%(q+1)!=1],nrow=p,ncol=q,byrow=TRUE)
      colnames(gamma0)<-Enames
      rownames(gamma0)<-Gnames
  est.res<-list(a=a0,beta=beta0, gamma=gamma0,b=b0,alpha=alpha0,coef=c(alpha0,a0,b0))
  class(est.res) = 'GEInfo'
  return(est.res)
}
