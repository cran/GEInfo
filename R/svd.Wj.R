svd.Wj <- function (Wj)
{
  n <- ifelse(is.vector(Wj), length(Wj), nrow(Wj))
  mat <- t(Wj) %*% Wj / n
  svd.res <- svd (mat)
  d.mat.sqrt.inverse <- diag(sqrt(svd.res$d) ^ -1, length(svd.res$d))
  inverse.mat <- svd.res$u %*% d.mat.sqrt.inverse
  Wj.tilde <- Wj %*% inverse.mat
  return (list(Wj.tilde = Wj.tilde, inverse.mat = inverse.mat))
}
