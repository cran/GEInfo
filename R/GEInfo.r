#' @title GEInfo approach with fixed tunings
#' @description Realize to estimate the GEInfo approach at fixed tunings.
#'  It is available for Linear, Logistic, and Poisson regressions.
#' @param E Observed matrix of E variables, of dimensions n x q.
#' @param G Observed matrix of G variables, of dimensions n x p.
#' @param Y Response variable, of length n. Quantitative for family="gaussian", or family="poisson" (non-negative counts). For family="binomial" should be a factor with two levels.
#' @param family Model type: one of ("gaussian", "binomial", "poisson").
#' @param kappa1 A user supplied kappa1.
#' @param kappa2 A user supplied kappa2.
#' @param lam1 A user supplied lambda1.
#' @param lam2 A user supplied lambda2.
#' @param tau A user supplied tau.
#' @param xi Tuning parameter of MCP penalty. Default is 6.
#' @param epsilon Tuning parameter of Ridge penalty which shrinks on the coefficients having prior information. Default is 0.
#' @param max.it Maximum number of iterations (total across entire path). Default is 500.
#' @param thresh Convergence threshold for group coordinate descent algorithm. The algorithm iterates until the change for each coefficient is less than thresh. Default is 1e-3.
#' @param S_G A user supplied vector, denoting the subscript of G variables which have prior information.
#' @param S_GE A user supplied matrix, denoting the subscript of G-E interactions which have prior information.
#'             The first and second columns of S_GE represent the subscript of G variable and the subscript of E variable, respectively.
#'             For example, S_GE = matrix( c(1, 2), ncol = 2), which indicates that the 1st G variable and the 2nd E variable have an interaction effect on Y.
#' @param Type_Y A vector of Type_Y prior information, having the same length with Y. Default is NULL.
#'              For family="gaussian", Type_Y is continuous. For family="binomial", Type_Y is binary.
#'              For family="poisson", Type_Y is a count vector.
#'              If users supply a Type_Y prior information, the function will use it to estimate a GEInfo model. If Type_Y=NULL,
#'              the function will incorporate the Type_S prior information S_G and S_GE to realize a GEInfo model.

#' @return An object of class "GEInfo" is returned, which is a list with the ingredients of the cross-validation fit.
#' \item{a}{Coefficient vector of length q for E variables.}
#' \item{b}{Coefficient vector of length (q+1)p for W (G variables and G-E interactions).}
#' \item{beta}{Coefficient vector of length p for G variables.}
#' \item{gamma}{Coefficient matrix of dimensions p*q for G-E interactions.}
#' \item{alpha}{Intercept.}
#' \item{coef}{A coefficient vector of length (q+1)*(p+1), including the estimates for \eqn{\alpha} (intercept), \eqn{a} (coefficients for all E variables), and \eqn{b} (coefficients for all G variables and G-E interactions).}

#' @details The function contains five tuning parameters, namely kappa1, kappa2, lambda1, lambda2, and tau.
#'  kappa1 and kappa2 are used to estimate model and select variables.
#'  lambda1 and lambda2 are used to calculate the prior-predicted response based on S_G and S_GE.
#'  tau is used for balancing between the observed response Y and the prior-predicted response.
#' @references Wang X, Xu Y, and Ma S. (2019). Identifying gene-environment interactions incorporating prior information. Statistics in medicine, 38(9): 1620-1633. \doi{10.1002/sim.8064}
#' @export
#' @examples
#' n <- 30; p <- 4; q <- 2
#' E <- MASS::mvrnorm(n, rep(0,q), diag(q))
#' G <- MASS::mvrnorm(n, rep(0,p), diag(p))
#' W <- matW(E, G)
#' alpha <- 0; a <- seq(0.4, 0.6, length=q);
#' beta <- c(seq(0.2, 0.5, length=2), rep(0, p-2))
#' vector.gamma <- c(0.8, 0.9, 0, 0)
#' gamma <- matrix(c(vector.gamma, rep(0, p*q - length(vector.gamma))), nrow=p, byrow=TRUE)
#' mat.b.gamma <- cbind(beta, gamma)
#' b <- as.vector(t(mat.b.gamma))              # coefficients of G and GE
#' Y <- alpha + E %*% a + W %*% b + rnorm (n, 0, 0.5)
#' S_G <- c(1)
#' S_GE <- cbind(c(1), c(1))
#' fit3 <- GEInfo(E, G, Y, family='gaussian', S_G=S_G,
#' S_GE=S_GE,kappa1 = 0.2,kappa2=0.2,lam1=0.2,lam2=0.2,tau=0.5)

GEInfo <- function(E, G, Y, family, S_G, S_GE, kappa1, kappa2, lam1, lam2, tau,
                   xi=6, epsilon=0, max.it=500, thresh=1e-03,
                   Type_Y=NULL) {
  q <- ncol(E)
  p <- ncol(G)
  Y.Info0 <- Type_Y
  ntau = length(tau)
  result.mat <- matrix(0, ncol = ntau, nrow = (p + 1) * (q + 1))
  critvec <- vector()
  N <- nrow(E)
  selec.kappa <- matrix(NA, nrow = ntau, ncol = 2)
  W <- matW(E, G)
  combi.X <- cbind(1, E, W)

  if (is.null(Type_Y)) {
    fit_CGEInfo <- CGEInfo(E, G, Y,family=family, xi=xi, epsilon=epsilon, max.it=max.it, thresh=thresh,
                           lam1=lam1, lam2=lam2, S_G=S_G, S_GE=S_GE)
    Type_Y <-stats::predict (object=fit_CGEInfo, Enew=E, Gnew=G, family=family,lam1=lam1,lam2=lam2)
  }

  Y.mod <- (1 - tau) * Y + tau * Type_Y
  fit.mod <- CGEInfo (E=E, G=G, Y=Y.mod, family = family,
                      xi=xi, epsilon=epsilon, max.it=max.it, thresh=thresh,
                      lam1=kappa1, lam2=kappa2)
  if(is.null(colnames(E))) {Enames<-paste0('E',1:q)} else { Enames<-colnames(E)}
  if(is.null(colnames(G))) {Gnames<-paste0('G',1:p)} else { Gnames<-colnames(G)}
   a0=fit.mod $ a
  names(a0)<-Enames

  b0=fit.mod $ b
  index.beta0<- seq(length(b0))%%(q+1)==1

  names(b0)[index.beta0]<-Gnames
  names(b0)[!index.beta0]<-paste0(rep(Gnames,each=q),'-',rep(Enames,p))
  beta0<-b0[index.beta0]

  alpha0=fit.mod $ alpha
  names(alpha0)<-'Intercept'
  gamma0<-fit.mod $ gamma
  colnames(gamma0)<-Enames
  rownames(gamma0)<-Gnames

  res<-list(a=a0,beta=beta0, gamma=gamma0,b=b0,alpha=alpha0,coef=c(alpha0,a0,b0))
  class(res) = 'GEInfo'
  return(res)
}
