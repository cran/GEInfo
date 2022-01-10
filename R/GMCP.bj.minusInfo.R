GMCP.bj.minusInfo <-
  function (Wj.Pen.tilde,
            bj.Pen,
            r.minus.j.Pen,
            kappa1,
            kappa2,
            xi) {
    n <-
      ifelse(is.vector(Wj.Pen.tilde),
             length(Wj.Pen.tilde),
             nrow(Wj.Pen.tilde))
    muj <- (t (Wj.Pen.tilde) %*% r.minus.j.Pen) / n
    norm.bj <- norm(bj.Pen, '2')
    df <- length(bj.Pen)
    if (kappa1 != 0) {
      ghat.part2 <-
        ifelse(norm.bj >= xi * sqrt(df) * kappa1,
               0,
               sqrt(df) * kappa1 - norm.bj / xi)
      ghat <- ifelse(norm.bj == 0, 1e+20, 1 + (1 / norm.bj) * ghat.part2)
    } else {
      ghat <- 1
    }
    vj <- ifelse (abs(muj) >= xi * kappa2 * ghat,
                  muj,
                  sapply (muj, Sfun1, kappa2) / (1 - 1 / (xi * ghat)))
    norm.vj <- norm (vj, '2')
    if (norm.vj >= xi * sqrt(df) * kappa1) {
      bjnew <- vj
    } else{
      bjnew <- (xi / (xi - 1)) * Sfun(vj, sqrt(df) * kappa1)
    }
    return(bj.Pen = bjnew)
  }
