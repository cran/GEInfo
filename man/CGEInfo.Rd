% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CGEInfo.R
\name{CGEInfo}
\alias{CGEInfo}
\title{CGEInfo and GEsgMCP approaches with fixed tunings}
\usage{
CGEInfo(
  E,
  G,
  Y,
  family,
  lam1,
  lam2,
  xi = 6,
  epsilon = 0,
  max.it = 500,
  thresh = 0.001,
  S_G = NULL,
  S_GE = NULL
)
}
\arguments{
\item{E}{Observed matrix of E variables, of dimensions n x q.}

\item{G}{Observed matrix of G variables, of dimensions n x p.}

\item{Y}{Response variable of length n. Quantitative for family="gaussian", or family="poisson" (non-negative count). For family="binomial" should be a factor with two levels.}

\item{family}{Model type: one of ("gaussian", "binomial", "poisson").}

\item{lam1}{A user supplied lambda1.}

\item{lam2}{A user supplied lambda2.}

\item{xi}{Tuning parameter of MCP penalty. Default is 6.}

\item{epsilon}{Tuning parameter of Ridge penalty which shrinks the coefficients of variables having prior information. Default is 0.}

\item{max.it}{Maximum number of iterations (total across entire path). Default is 500.}

\item{thresh}{Convergence threshold for group coordinate descent algorithm. The algorithm iterates until the change for each coefficient is less than thresh. Default is 1e-3.}

\item{S_G}{A user supplied vector, denoting the subscript of G variables which have prior information. Default is NULL.}

\item{S_GE}{A user supplied matrix, denoting the subscript of G-E interactions which have prior information.
The first and second columns of S_GE represent the subscript of G variable and the subscript of E variable, respectively.
For example, S_GE = matrix( c(1, 2), ncol = 2), which indicates that the 1st G and the 2nd E variables have an interaction effect on Y.
Default is NULL. If both S_G and S_GE are NULL, no prior information is incorporated in the model, in which case function CGEInfo realizes GEsgMCP approach.}
}
\value{
An object of class "GEInfo" is returned, which is a list including the estimation results at fixed tunings.
\item{a}{Coefficient vector of length q for E variables.}
\item{b}{Coefficient vector of length (q+1)p for W (G variables and G-E interactions).}
\item{beta}{Coefficient vector of length p for G variables.}
\item{gamma}{Coefficient matrix of dimensions p*q for G-E interactions.}
\item{alpha}{Intercept.}
\item{coef}{A coefficient vector of length (q+1)*(p+1), including the estimates for \eqn{\alpha} (intercept), \eqn{a} (coefficients for all E variables), and \eqn{b} (coefficients for all G variables and G-E interactions).}
}
\description{
Realize to estimate CGEInfo and GEsgMCP approaches at fixed tunings.
}
\examples{
n <- 30; p <- 5; q <- 2
E <- MASS::mvrnorm(n, rep(0,q), diag(q))
G <- MASS::mvrnorm(n, rep(0,p), diag(p))
W <- matW(E, G)
alpha <- 0; a <- seq(0.4, 0.6, length=q);
beta <- c(seq(0.2, 0.5, length=3),rep(0, p-3))  # coefficients of G variables
vector.gamma <- c(0.8, 0.5, 0, 0)
gamma <- matrix(c(vector.gamma, rep(0, p*q - length(vector.gamma))), nrow=p, byrow=TRUE)
mat.b.gamma <- cbind(beta, gamma)
b <- as.vector (t(mat.b.gamma))              # coefficients of G and G-E interactions
Y <- alpha + E \%*\% a + W \%*\% b + rnorm (n, 0, 0.5)
S_G <- c(1)
S_GE <- cbind(c(1), c(1))
fit1 <- CGEInfo(E, G, Y,family='gaussian', S_G=S_G, S_GE=S_GE,lam1=0.4,lam2=0.4)
}
\references{
Wang X, Xu Y, and Ma S. (2019). Identifying gene-environment interactions incorporating prior information. Statistics in medicine, 38(9): 1620-1633. \doi{10.1002/sim.8064}
}
