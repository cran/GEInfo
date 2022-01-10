#' @title Calculate matrix W
#' @description Calculate observed matrix W for all G variables and G-E interactions. Denote Wj as the n x (q+1) sub-matrix of W corresponding the jth G variable.
#'  The first column of Wj is the observation vector of the jth G variable, and the rest q columns of Wj are observations of G-E interactions.
#' @param E Observed matrix of E variables, of dimension n x q.
#' @param G Observed matrix of G variables, of dimensions n x p.
#' @return A matrix of dimension n x [p(q+1)].
#' @export
#' @examples
#' n <- 30; q <- 3; p <- 5;
#' E <- MASS::mvrnorm (n, rep (0, q), diag (q))
#' G <- MASS::mvrnorm (n, rep (0, p), diag (p))
#' W <- matW (E, G)

matW <- function(E, G)
{
  Estar <- cbind(1, E)
  p <- ncol(G)
  q <- ncol(E)
  n <- nrow(E)
  W <- matrix(0, n, p * (q + 1))
  for (j in 1:p) {
    ind <- ((q + 1) * (j - 1) + 1):((q + 1) * j)
    Gj <- G[, j]
    Wj <- Gj * Estar
    W[, ind] <- Wj
  }
  return(W = W)
}

