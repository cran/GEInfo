update.para.svd.c <-
  function (Eols,
            Wols,
            Yols,
            AS_c,
            alpha0,
            a0,
            b0,
            kappa1,
            kappa2,
            xi,
            epsilon,
            S_G,
            S_GE) {
    q<-ncol(Eols)-1 # Eols includes intercept and E
    p<-ncol(Wols)/(q+1)
    r <- Yols - Eols %*% c(alpha0, a0) - Wols %*% b0
    groupid <- ceiling(AS_c / (q + 1))
    inactive.G <- unique(groupid)
    for (j in inactive.G) {
      jind <- AS_c[groupid == j]
      if (min(jind) %% (q + 1) == 1) {
        Wj <- Wols[, jind]
        svd.res <- svd.Wj(Wj)
        Wj.tilde <- svd.res$Wj.tilde
        inverse.mat <-
          svd.res$inverse.mat
        bj <- b0[jind]
        r.minus.j <-  r + crossprod(t(Wj), bj)
        fit <- GMCP.bj(Wj.tilde, bj, r.minus.j, kappa1, kappa2, xi)
        bj.new <- inverse.mat %*% fit
        b0[jind] <- bj.new
        r <- r.minus.j - crossprod(t(Wj), bj.new)
      } else  {
        Wj <- Wols[, jind]
        svd.res <- svd.Wj(Wj)
        Wj.tilde <- svd.res$Wj.tilde
        inverse.mat <-
          svd.res$inverse.mat
        bj <- b0[jind]
        r.minus.j <-  r +  crossprod(t(Wj), bj)
        fit <-
          GMCP.bj.minusInfo (Wj.tilde, bj, r.minus.j, kappa1, kappa2, xi)
        bj.new <- inverse.mat %*% fit
        b0[jind] = bj.new
        r <- r.minus.j - crossprod(t(Wj), bj.new)
      }
    }
    return(b = round(b0, 4))
  }
