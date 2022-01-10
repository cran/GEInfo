update.para.svd  <-
  function (Eols,
            Wols,
            Yols,
            active.set,
            alpha0,
            a0,
            b0,
            kappa1,
            kappa2,
            xi,
            epsilon,
            S_G,
            S_GE) {
    q <- ncol(Eols)-1  # Eols includes intercept and E
    p <- ncol(Wols)/(q+1)
    index <- active.set
    r <- Yols - crossprod(t(Wols[, index]), b0[index])
    a0 <- update.a (Eols, r)$a0
    r <- r - Eols %*% a0
    groupid <- ceiling(active.set / (q + 1))
    active.G <- unique(groupid)
    if (is.null (S_G) & is.null (S_GE)) {
      for (j in active.G) {
        jind0 <- active.set[groupid == j]
        jind <- union(j * (q + 1) - q, jind0)
        if (length(jind) == 1) {
          Wj1 <- Wols[, jind]
          bj1 <- b0[jind]
          r.minus.j1 <- r + Wj1 * bj1
          muj1 <- Wj1 %*% r.minus.j1 / length(Yols)
          wt <- sum(Wj1 ^ 2) / length(Yols)
          fit <-
            ifelse (abs(muj1) >= xi * kappa1 * wt,
                    muj1 / wt,
                    Sfun1(muj1, kappa1) / (wt - 1 / xi))
          b0[jind] <- fit
          r <- r.minus.j1 - Wj1 * b0[jind]
        } else if (length(jind) > 1) {
          Wj <- Wols[, jind]
          svd.res <- svd.Wj(Wj)
          Wj.tilde <- svd.res$Wj.tilde
          inverse.mat <- svd.res$inverse.mat
          bj <- b0[jind]
          r.minus.j <-  r + Wj %*% bj
          fit <- GMCP.bj(Wj.tilde, bj, r.minus.j, kappa1, kappa2, xi)
          bj.new <- inverse.mat %*% fit
          b0[jind] <- bj.new
          r <- r.minus.j - Wj %*% bj.new
        }
      }
    } else {
      S_H <- unique(S_GE[, 1])
      S_G_all <- sort(union(S_G, S_H))
      if (any(S_GE[, 2] > q)) {
        message ('ERROR: The elements of the 2nd colume of S_GE should not exceed q')
      }   else if (any(S_G > p)) {
        message ('ERROR: The elements of S_G should not exceed p')
      }   else if (any(S_GE[1, ] > p)) {
        message ('ERROR: The elements of the 1st colume of S_GE should not exceed p')
      } else {
        beta.ind <- (q + 1) * (S_G_all - 1) + 1
        jind <- S_GE[, 1]
        kind <- S_GE[, 2]
        gamma.ind <- (q + 1) * (jind - 1) + 1 + kind
        MLE.ind <- sort(c(beta.ind, gamma.ind))
        Wj.info <- Wols[, MLE.ind]
        bj.MLE.old <- b0[MLE.ind]
        r.minus.MLE <- r + crossprod(t(Wj.info), bj.MLE.old)
        if (length(MLE.ind) == 1) {
          b0[MLE.ind] <- stats::coef(stats::lm(r.minus.MLE ~ Wj.info))[-1]
        } else {
          fit <- MASS::lm.ridge(r.minus.MLE ~ Wj.info, lambda = epsilon)
          b0[MLE.ind] = stats::coef(fit)[-1]
        }
        r <- r.minus.MLE - crossprod(t(Wj.info), b0[MLE.ind])
      }
      set.pen <- setdiff (active.set, MLE.ind)
      groupid0 <- ceiling(set.pen / (q + 1))
      active.G0 <- unique(groupid0)
      for (j in active.G0) {
        jind <- set.pen[groupid0 == j]
        if (length(jind) == 1) {
          if (jind %% (q + 1) ==  1) {
            Wj1 <- Wols[, jind]
            bj1 <- b0[jind]
            r.minus.j1 <- r + crossprod(t(Wj1), bj1)
            muj1 <- t(Wj1) %*% r.minus.j1 / length(Yols)
            wt <- sum(Wj1 ^ 2) / length(Yols)
            fit <-
              ifelse (abs(muj1) >= xi * kappa1 * wt,
                      muj1 / wt,
                      Sfun1(muj1, kappa1) / (wt - 1 / xi))
            b0[jind] <- fit
            r <- r.minus.j1 - Wj1 * b0[jind]
          }
        }
        if (min(jind) %% (q + 1) != 1) {
          Wj <- Wols[, jind]
          svd.res <- svd.Wj(Wj)
          Wj.tilde <- svd.res$Wj.tilde
          inverse.mat <- svd.res$inverse.mat
          bj <- b0[jind]
          r.minus.j <-  r +  crossprod(t(Wj), bj)
          fit <-
            GMCP.bj.minusInfo (Wj.tilde, bj, r.minus.j, kappa1, kappa2, xi)
          bj.new <- inverse.mat %*% fit
          b0[jind] = bj.new
          r <- r.minus.j - crossprod(t(Wj), bj.new)
        }
        if (min(jind) %% (q + 1) == 1 & length(jind) > 1)  {
          Wj <- Wols[, jind]
          svd.res <- svd.Wj(Wj)
          Wj.tilde <- svd.res$Wj.tilde
          inverse.mat <- svd.res$inverse.mat
          bj <- b0[jind]
          r.minus.j <-  r + Wj %*% bj
          fit <- GMCP.bj(Wj.tilde, bj, r.minus.j, kappa1, kappa2, xi)
          bj.new <- inverse.mat %*% fit
          b0[jind] <- bj.new
          r <- r.minus.j - Wj %*% bj.new
        }
      }
    }
    res <- list(a = round(a0, 4), b = round(b0, 4))
    return(res)
  }
