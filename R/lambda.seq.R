lambda.seq <- function(E,
                       G,
                       Y,
                       family,
                       lam1 = NULL,
                       lam2 = NULL) {
  n <- length(Y)
  W <- matW(E, G)
  p <- ncol(G)
  q <- ncol(E)
  if (is.null(lam1))
  {
    if (family == 'poisson') {
      lam1 <- c(0.01, 0.05, 0.1, 0.5, 0.8, 1, 5, 7, 10, 50)
    }
    if (family == 'gaussian')
    {
      fit <- stats::lm(Y ~ E)
      residual.a <- Y - fit$fitted.values
      value1 = 0
      for (j in seq(p))
      {
        ind <- (j - 1) * (q + 1) + 1:(q + 1)
        value1[j] <-
          norm(t(W[, ind]) %*% residual.a , '2') / (n * sqrt (q + 1))
      }
      lam1.max <-  max(value1)
      lam1.min <- 0.006 * lam1.max
      lam1 <- 10 ^ seq(log10(lam1.min), log10(lam1.max), length = 10)
    }
    if (family == 'binomial')
    {
      fit <- stats::glm(ifelse(Y >= 0.5, 1, 0) ~ E, family = 'binomial')
      residual.a <- Y - fit$fitted.values
      value1 = 0
      for (j in seq(p))
      {
        ind <- (j - 1) * (q + 1) + 1:(q + 1)
        value1[j] <-
          norm(t(W[, ind]) %*% residual.a , '2') / (n * sqrt (q + 1))
      }
      lam1.max <-  max(value1) * 3
      lam1.min <- 0.08 * lam1.max
      lam1 <- 10 ^ seq(log10(lam1.min), log10(lam1.max), length = 10)
    }
  }

  if (is.null(lam2))
  {
    if (family == 'poisson')
    {
      lam2 <- c(0.01, 0.05, 0.1, 0.5, 0.8, 1, 5, 7, 10, 50)
    }
    if (family == 'gaussian')
    {
      fit <- stats::lm(Y ~ E)
      residual.a <- Y - fit$fitted.values
      lam2.max <- max(abs (t(W) %*% residual.a / n))
      lam2.min <- 0.006 * lam2.max
      lam2 <- 10 ^ seq(log10(lam2.min), log10(lam2.max), length = 10)
    }
    if (family == 'binomial')
    {
      fit <- stats::glm(ifelse(Y >= 0.5, 1, 0) ~ E, family = 'binomial')
      residual.a <- Y - fit$fitted.values
      lam2.max <- max(abs (t(W) %*% residual.a / n)) * 3
      lam2.min <- 0.06 * lam2.max
      lam2 <- 10 ^ seq(log10(lam2.min), log10(lam2.max), length = 10)
    }
  }
  return(list(lambda1 = lam1, lambda2 = lam2))
}
