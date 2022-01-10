GMCP.bj <- function (Wj.tilde, bj, r.minus.j, kappa1, kappa2, xi) {
  n <- nrow(Wj.tilde)
  df <- ncol(Wj.tilde)
  muj <- (t (Wj.tilde) %*% r.minus.j) / n
  norm.bj <- norm(bj, '2')
  if (kappa1 != 0) {
    ghat.part2 <-
      ifelse(norm.bj >= xi * sqrt(df) * kappa1,
             0,
             sqrt(df) * kappa1 - norm.bj / xi)
    ghat <- ifelse(norm.bj == 0, 1e+20, 1 + (1 / norm.bj) * ghat.part2)
  } else {
    ghat <- 1
  }
  vj <- muj
  vj[-1] <- ifelse (abs(muj[-1]) >= xi * kappa2 * ghat,
                    muj[-1],
                    sapply (muj[-1], Sfun1, kappa2) / (1 - 1 / (xi * ghat)))
  norm.vj <- norm (vj, '2')
  if (norm.vj >= xi * sqrt(df) * kappa1) {
    bjnew <- vj
  } else{
    bjnew <- (xi / (xi - 1)) * Sfun(vj, sqrt(df) * kappa1)
  }
  return(bj = bjnew)
}
