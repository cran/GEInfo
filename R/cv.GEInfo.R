#' @title  Cross-validation for GEInfo
#' @description Does k-fold cross-validation for GEInfo approach,
#'  which adaptively accommodates the quality of the prior information and automatically detects the false information.
#'  Tuning parameters are chosen based on a user given criterion.
#' @param E Observed matrix of E variables, of dimensions n x q.
#' @param G Observed matrix of G variables, of dimensions n x p.
#' @param Y Response variable, of length n. Quantitative for family="gaussian", or family="poisson" (non-negative counts). For family="binomial" should be a factor with two levels.
#' @param family Model type: one of ("gaussian", "binomial", "poisson").
#' @param nfolds Number of folds. Default is 3.
#'              Although nfolds can be as large as the sample size n (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3
#' @param xi Tuning parameter of MCP penalty. Default is 6.
#' @param epsilon Tuning parameter of Ridge penalty which shrinks on the coefficients having prior information. Default is 0.
#' @param max.it Maximum number of iterations (total across entire path). Default is 500.
#' @param thresh Convergence threshold for group coordinate descent algorithm. The algorithm iterates until the change for each coefficient is less than thresh. Default is 1e-3.
#' @param criterion Criterion used for tuning selection via cross-validation. Currently five options: MSE, AIC, BIC, EBIC, GCV. Default is BIC. See Details.
#' @param S_G A user supplied vector, denoting the subscript of G variables which have prior information.
#' @param S_GE A user supplied matrix, denoting the subscript of GE interactions which have prior information.
#'             The first and second columns of S_GE represent the subscript of G variable and the subscript of E variable, respectively.
#'             For example, S_GE = matrix( c(1, 2), ncol = 2), which indicates that the 1st G variable and the 2nd E variables have an interaction effect on Y.
#' @param Type_Y A vector of Type_Y prior information, having the same length with Y. Default is NULL.
#'              For family="gaussian", Type_Y is continuous. For family="binomial", Type_Y is binary.
#'              For family="poisson", Type_Y is count.
#'              If users supply a Type_Y prior information, this function will use it to estimate a GEInfo model. If Type_Y=NULL,
#'              the function will incorporate the prior information included in S_G and S_GE to realize a GEInfo model.

#' @param kappa1 A user supplied kappa1 sequence. Default is kappa1=NULL.
#'             Typical usage is to have the program compute its own kappa1 sequence. Supplying a value of kappa1 overrides this. See Details.
#' @param kappa2 A user supplied kappa2 sequence. Default is kappa2=NULL.
#'             Typical usage is to have the program compute its own kappa2 sequence. Supplying a value of kappa2 overrides this. See Details.

#' @param lam1 A user supplied lambda1 sequence. Default is lam1=NULL.
#'             Typical usage is to have the program compute its own lambda1 sequence. Supplying a value of lam1 overrides this. See Details.
#' @param lam2 A user supplied lambda2 sequence. Default is lam2=NULL.
#'             Typical usage is to have the program compute its own lambda1 sequence. Supplying a value of lam2 overrides this. See Details.
#' @param tau A user supplied tau sequence ranging from 0 to 1.
#'            Default is tau = c (0, 0.25,0.5,0.75,1). See Details.
#' @details The function contains five tuning parameters, namely kappa1, kappa2, lambda1, lambda2, and tau.
#'  kappa1 and kappa2 are used to estimate model and select variables.
#'  lambda1 and lambda2 are used to calculate the prior-predicted response based on S_G and S_GE.
#'  tau is used for balancing between the observed response Y and the prior-predicted response.
#'  When tau=0 and tau=1, this function realizes cross-validation for GEsgMCP and CGEInfo approaches, respectively.
#'
#' In order to select the optimal tuning combination, there are five criteria available, which are MSE, AIC, BIC, GCV, and EBIC. Let L be the loss function of the model,
#' \eqn{MSE=L}, \eqn{AIC=2L+2df}, \eqn{BIC=2L+ln(n)df}, \eqn{GCV=2L/(1-df/n)^2},
#' and \eqn{EBIC=2L+ln(n)df + 2df ln(nvar) (1-ln(n)/(2ln(nvar)))}.
#'  In most cases, BIC is a good choice. In the case of high dimension, EBIC criterion is recommended first,
#'  which has demonstrated satisfactory performance in high-dimensional studies.

#' @return An object of class "GEInfo" is returned, which is a list with the ingredients of the cross-validation fit.
#'\item{coef.all.tau}{A matrix of coefficients, of dimensions (p+1)(q+1) x length(tau).}
#'\item{best.tuning}{A list containing the optimal tau, kappa1, and kappa2.}
#' \item{a}{Coefficient vector of length q for E variables.}
#' \item{beta}{Coefficient vector of length p for E variables.}
#' \item{gamma}{Coefficient matrix of dimensions p*q for G-E interactions.}
#' \item{b}{Coefficient vector of length (q+1)p for W (G variables and G-E interactions).}
#' \item{alpha}{Intercept.}
#' \item{coef}{A coefficient vector of length (q+1)(p+1), including the estimates for \eqn{\alpha} (intercept), \eqn{a} (coefficients for all E variables), and \eqn{b} (coefficients for all G variables and G-E interactions).}
#' \item{nvar}{Number of non-zero coefficients at the best tunings.}
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
#' b <- as.vector(t(mat.b.gamma))
#' Y <- alpha + E %*% a + W %*% b + rnorm (n, 0, 0.5)
#' S_G <- c(1)
#' S_GE <- cbind(c(1), c(1))
#' fit4 <- cv.GEInfo(E, G, Y, family='gaussian', S_G=S_G,
#'  S_GE=S_GE,lam1=0.4,lam2=0.4,kappa1 = 0.4,kappa2=0.4,tau=0.5)

cv.GEInfo <- function(E, G, Y, family, S_G, S_GE,
                      nfolds=3, xi=6, epsilon=0, max.it=500, thresh=1e-03,
                      criterion='BIC', Type_Y=NULL, kappa1=NULL, kappa2=NULL, lam1=NULL, lam2=NULL, tau=c(0, 0.25,0.5,0.75,1)) {
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
    fit_CGEInfo <- cv.CGEInfo(E, G, Y,family=family, nfolds=nfolds, xi=xi, epsilon=epsilon, max.it=max.it, thresh=thresh,
                              criterion=criterion, lam1=lam1, lam2=lam2, S_G=S_G, S_GE=S_GE)
    Type_Y <-stats::predict (object=fit_CGEInfo, Enew=E, Gnew=G, family=family) }
  for (i in seq(ntau)) {
    if(tau[i]==1 & is.null(Y.Info0)) {
      fit.mod <-fit_CGEInfo
    } else {
      Y.mod <- (1 - tau[i]) * Y + tau[i] * Type_Y
      fit.mod <- cv.CGEInfo (E=E, G=G, Y=Y.mod, family = family,
                             nfolds=nfolds, xi=xi, epsilon=epsilon, max.it=max.it, thresh=thresh,
                             criterion=criterion, lam1=kappa1, lam2=kappa2)
    }
    selec.kappa[i, ] <- fit.mod$best.tuning
    result.mat[, i] <- stats::coef(fit.mod)
    index.nonzero <- which(stats::coef(fit.mod) != 0)
    eta.pred <- crossprod(t(combi.X[, index.nonzero]), stats::coef(fit.mod)[index.nonzero])
    df <- sum(stats::coef(fit.mod)[-1] != 0)
    critvec[i] <- criteria.pred(test_Y = Y, eta.pred = eta.pred,
                                p=p, q=q, family=family, criterion=criterion, df=df) / N
  }
  min.ind <- max(which(critvec == min(critvec), arr.ind = TRUE))
  best.tau <- tau[min.ind]
  best.kappa <- selec.kappa [min.ind, ]
  if(is.null(colnames(E))) {Enames<-paste0('E',1:q)} else { Enames<-colnames(E)}
  if(is.null(colnames(G))) {Gnames<-paste0('G',1:p)} else { Gnames<-colnames(G)}

  best.est <- round(result.mat[, min.ind], 4)
  a.est<- best.est[2:(q+1)]
  names(a.est)<-Enames

  b.est<-best.est[-(1:(q+1))]
  index.beta0<- seq(length(b.est))%%(q+1)==1
  names(b.est)[index.beta0]<-Gnames
  names(b.est)[!index.beta0]<-paste0(rep(Gnames,each=q),'-',rep(Enames,p))
  beta.est <- b.est [index.beta0]
  gamma.est <- matrix(b.est[!index.beta0],ncol=q,nrow=p,byrow = T)
  colnames(gamma.est)<-Enames
  rownames(gamma.est)<-Gnames
  intercept<-best.est[1]
  names(intercept)<-'Intercept'
  coef0<-c(intercept,a.est,b.est)
  res<-list(coef.all.tau=result.mat, best.tuning=list(tau=best.tau, kappa1=best.kappa[1], kappa2=best.kappa[2]),
            alpha=intercept, a=a.est, beta=beta.est,
            gamma=gamma.est, b=b.est, coef=coef0, nvar=sum(coef0[-1]!=0))
  class(res) = 'GEInfo'
  return(res)
}
