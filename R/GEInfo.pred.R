GEInfo.pred <-
  function (obj,
            test_E,
            test_G,
            test_Y,
            nlambda1,
            nlambda2,
            family,
            criterion) {
    n.test <- nrow(test_E)
    q <- ncol(test_E)
    p <- ncol(test_G)
    test_W <- matW(test_E, test_G)
    a.est <- obj$a
    b.est <- obj$b
    alpha.est <- obj$alpha
    df.mat <- matrix(0, nlambda1, nlambda2)
    mat.pred <- matrix(0, nlambda1, nlambda2)
    for (i in 1:nlambda1) {
      for (j in 1:nlambda2) {
        df.col <- sum(obj$est.all[i, j, -1] != 0)
        df.mat[i, j] <- df.col
      }
    }
    combi.X <- cbind(1, test_E, test_W)
    for (i in seq(nlambda1))
    {
      for (j in seq(nlambda2))
      {
        eta.ij <- combi.X %*% obj$est.all [i, j, ]
        mat.pred[i, j] <-
          criteria.pred(
            eta.pred = eta.ij,
            test_Y = test_Y,
            p = p,
            q = q,
            family = family,
            criterion = criterion,
            df = df.mat[i, j]
          )
      }
    }
    return(mat.pred)
  }
