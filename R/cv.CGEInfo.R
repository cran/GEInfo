#' @title Cross-validation for CGEInfo
#' @description Does k-fold cross-validation for CGEInfo, returns the estimation results at best tunings, and produces a heatmap for the identification results.
#' @param E Observed matrix of E variables, of dimensions n x q.
#' @param G Observed matrix of G variables, of dimensions n x p.
#' @param Y Response variable, of length n. Quantitative for family="gaussian", or family="poisson" (non-negative counts). For family="binomial" should be a factor with two levels.
#' @param family Model type: one of ("gaussian", "binomial", "poisson").
#' @param nfolds Number of folds. Default is 3.
#' Although nfolds can be as large as the sample size n (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3.
#' See Details.
#' @param xi Tuning parameter of MCP penalty. Default is 6.
#' @param epsilon Tuning parameter of Ridge penalty which shrinks the coefficients having prior information. Default is 0.
#' @param max.it Maximum number of iterations (total across entire path). Default is 500.
#' @param thresh Convergence threshold for group coordinate descent algorithm. The algorithm iterates until the change for each coefficient is less than thresh. Default is 1e-3.
#' @param criterion Criterion used for cross-validation. Currently five options: MSE, AIC, BIC, EBIC, GCV. Default is BIC. See Details.
#' @param lam1 A user supplied  lambda1 sequence.
#'            Typical usage is to have the program compute its own lambda1 sequence. Supplying a value of lam1 overrides this. Default is lam1=NULL.
#' @param lam2 A user supplied  lambda2 sequence. Default is lam2=NULL.
#'            Typical usage is to have the program compute its own lambda2 sequence. Supplying a value of lam2 overrides this. Default is lam2=NULL.
#' @param S_G A user supplied vector, denoting the subscript of G variables which have prior information. Default is NULL. See Details.
#' @param S_GE A user supplied matrix, denoting the subscript of G-E interactions which have prior information.
#'             The first and second columns of S_GE represent the subscript of G variable and the subscript of E variable, respectively.
#'             For example, S_GE = matrix( c(1, 2), ncol = 2), indicating that the 1st G variable and the 2nd E variables have an interaction effect on Y.
#'             Default is NULL. If both S_G and S_GE are NULL, no prior information is incorporated in the model, in which case this function realizes GEsgMCP approach. See Details.
#'
#' @return An object of class "GEInfo" is returned, which is a list with the ingredients of the cross-validation fit.
#' \item{best.tuning}{A vector of length 2, containing the best lambda1 and lambda2 selected by cross-validation.}
#' \item{a}{Coefficient vector of length q for all E variables.}
#' \item{beta}{Coefficient vector of length p for all G variables.}
#' \item{gamma}{Coefficient matrix of dimensions p*q for G-E interactions.}
#' \item{b}{Coefficient vector of length (q+1)*p for W (G variables and G-E interactions).}
#' \item{alpha}{Intercept.}
#' \item{coef}{A coefficient vector of length (q+1)*(p+1), including the estimates for \eqn{\alpha} (intercept), \eqn{a} (coefficients for all E variables), and \eqn{b} (coefficients for all G variables and G-E interactions).}
#' \item{nvar}{Number of non-zero coefficients at the best tunings.}

#' @details The function calls CGEInfo nfolds times, each time leaving out 1/nfolds of the data. The cross-validation error is based on the user given "criterion".
#' cv.CGEInfo supports to construct two methods: GEInfo and GEsgMCP, depending on whether S_G and S_GE are NULL.
#' When either S_G or S_GE is not NULL, CGEInfo approach is realized, which completely trusts the prior information.
#' Otherwise, GEsgMCP approach is constructed, in which no prior information is incorporated.
#'
#' In order to select the optimal tunings, there are five criteria available, which are MSE, AIC, BIC, GCV, and EBIC. Let L be the loss function of the model,
#' \eqn{MSE=L}, \eqn{AIC=2L+2df}, \eqn{BIC=2L+ln(n)df}, \eqn{GCV=2L/(1-df/n)^2},
#' and \eqn{EBIC=2L+ln(n)df + 2df ln(nvar) (1-ln(n)/(2ln(nvar)))}.
#'  In most cases, BIC is a good choice. In the case of high dimension, EBIC criterion is recommended first,
#'  which has demonstrated satisfactory performance in high-dimensional studies.
#' @references Wang X, Xu Y, and Ma S. (2019). Identifying gene-environment interactions incorporating prior information. Statistics in medicine, 38(9): 1620-1633. \doi{10.1002/sim.8064}
#' @export
#' @examples
#' n <- 30; p <- 5; q <- 2
#' E <- MASS::mvrnorm(n, rep(0,q), diag(q))
#' G <- MASS::mvrnorm(n, rep(0,p), diag(p))
#' W <- matW(E, G)
#' alpha <- 0; a <- seq(0.4, 0.6, length=q);
#' beta <- c(seq(0.2, 0.5, length=3),rep(0, p-3))
#' vector.gamma <- c(0.8, 0.5, 0, 0)
#' gamma <- matrix(c(vector.gamma, rep(0, p*q - length(vector.gamma))), nrow=p, byrow=TRUE)
#' mat.b.gamma <- cbind(beta, gamma)
#' b <- as.vector (t(mat.b.gamma))
#' Y <- alpha + E %*% a + W %*% b + rnorm (n, 0, 0.5)
#' S_G <- c(1)
#' S_GE <- cbind(c(1), c(1))
#' fit2 <- cv.CGEInfo(E, G, Y,family='gaussian', S_G=S_G, S_GE=S_GE,lam1=0.4,lam2=0.4)


cv.CGEInfo  <- function(E, G, Y, family, nfolds=3, xi=6, epsilon=0, max.it=500, thresh=1e-03,
                        criterion='BIC', lam1=NULL, lam2=NULL, S_G=NULL, S_GE=NULL) {
  n <- length(Y)
  p<-ncol(G)
  q<-ncol(E)
  lam <- lambda.seq(E, G, Y, family, lam1 = lam1, lam2 = lam2)
  lambda1 <- lam$lambda1
  lambda2 <- lam$lambda2
  nlambda1 <- length(lambda1)
  nlambda2 <- length(lambda2)
  predmat <- matrix(0, nrow = nlambda1, ncol = nlambda2)
  foldid = sample(rep(seq(nfolds), length = n))
  outlist <- list()
  for (k in seq (nfolds)) {
    which = foldid == k
    Y_sub <- Y[!which]
    G_sub <- G[!which,]
    E_sub <- E[!which,]
    if (family == 'poisson') {
      fit <- CGEInfo.poisson(E=E_sub, G=G_sub, Y=Y_sub, family='poisson', lam1=lambda1, lam2=lambda2,
                             xi=xi, epsilon=epsilon, max.it=max.it, thresh=thresh,S_G=S_G, S_GE=S_GE)
    } else {
      fit <- CGEInfo1 (E=E_sub, G=G_sub, Y=Y_sub, family=family, lam1=lambda1, lam2=lambda2,
                      xi=xi, epsilon=epsilon, max.it=max.it, thresh=thresh, S_G=S_G, S_GE=S_GE)
    }
    outlist[[k]] <- fit
    preds<-GEInfo.pred(obj=fit, test_E=E[which,], test_G=G[which,], test_Y=Y[which], nlambda1, nlambda2, family, criterion)
    preds0 <- unlist(preds)
    preds1 <- matrix(preds0, nrow = nlambda1, ncol = nlambda2)
    predmat <- predmat + preds1
  }
  critmat <- predmat / n
  matind <- which(critmat == min(critmat), arr.ind = TRUE)
  lam1ind <- max(matind[, 1])
  lam2ind <- max(matind[which(matind[, 1] == lam1ind), 2])
  best.tuning <- c(lambda1[lam1ind], lambda2[lam2ind])
  if(family=='poisson') {
    fit.final<-CGEInfo.poisson(E=E, G=G, Y=Y, family=family, lam1=lambda1[lam1ind], lam2=lambda2[lam2ind], xi=xi,
                               epsilon=epsilon, max.it=max.it, thresh=thresh, S_G=S_G, S_GE=S_GE)
  } else {
    fit.final<-CGEInfo1(E=E, G=G,Y=Y, family=family, lam1=lambda1[lam1ind], lam2=lambda2[lam2ind], xi=xi,
                       epsilon=epsilon, max.it=max.it, thresh=thresh, S_G=S_G, S_GE=S_GE)
  }
  if(is.null(colnames(E))) {Enames<-paste0('E',1:q)} else { Enames<-colnames(E)}
  if(is.null(colnames(G))) {Gnames<-paste0('G',1:p)} else { Gnames<-colnames(G)}

  a.est<-as.vector(unlist(fit.final$a))
  names(a.est)<-Enames

  b.est <- as.vector(unlist(fit.final$b))
  index.beta0<- seq(length(b.est))%%(q+1)==1
  names(b.est)[index.beta0]<-Gnames
  names(b.est)[!index.beta0]<-paste0(rep(Gnames,each=q),'-',rep(Enames,p))

  beta.est <- b.est [index.beta0]
  gamma.est <- matrix(b.est[!index.beta0],ncol=q,nrow=p,byrow = T)
  colnames(gamma.est)<-Enames
  rownames(gamma.est)<-Gnames

  intercept<-fit.final$alpha
  names(intercept)<-'Intercept'
  coef0<-c(intercept,a.est,b.est)
  res<-list(best.tuning=best.tuning,a=a.est,beta=beta.est,gamma=gamma.est,
            b=b.est, alpha=fit.final$alpha,
            coef=coef0, nvar=sum(coef0[-1]!=0))
  class(res) = 'GEInfo'
  return(res)
}
