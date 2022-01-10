criteria.pred <-
  function(test_Y,
           eta.pred,
           family,
           criterion,
           q,
           p,
           df = NULL) {
    N <- length(test_Y)
    num.cov <- p + q * p + q
    r.EBIC <- 1 - 1 / (2 * log(num.cov) / log(N))
    if (family == 'gaussian') {
      error.test <- test_Y - eta.pred
      loss.fun <- sum(error.test ^ 2) / 2
      return(switch(
        criterion,
        MSE = loss.fun,
        EBIC = 2 * loss.fun + df * log(N)  + 2 * df * r.EBIC * log(num.cov),
        AIC = 2 * loss.fun + 2 * df,
        BIC = 2 * loss.fun + log(N) * df,
        GCV = 2 * loss.fun / (1 - df / N) ^ 2
      ))
    }
    if (family == 'poisson') {
      lambda <- exp(eta.pred)
      pred_Y <- lambda
      logloss <- -(test_Y %*% eta.pred - sum(exp(eta.pred)))
      return(switch(
        criterion,
        MSE = logloss,
        EBIC = 2*logloss + df * log(N) /N + 2 * df * r.EBIC *
          log(num.cov)/N,
        AIC =  logloss + 2 * df/N,
        BIC =  logloss + log(N) * df/N,
        GCV =  logloss / (1 - df / N) ^ 2
      ))
    }
    if (family == 'binomial')
    {
      prob.pred <- 1 / (1 + exp(-eta.pred))
      pred_Y <- ifelse(prob.pred >= 0.5, 1, 0)
      cont <- table(test_Y, pred_Y)
      logloss <- -sum(ifelse(
        prob.pred == 1 | prob.pred == 0,
        -10 ^ 10,
        test_Y * log(prob.pred) + (1 - test_Y) * log(1 - prob.pred)
      ))
      return(switch(
        criterion,
        MSE =1- (cont[1,1]+cont[2,2]) / N,
        EBIC =  2*logloss + df * log(N)/N + 2 * df * r.EBIC *
          log(num.cov)/N,
        AIC = logloss + 2 * df/N,
        BIC =  logloss + log(N) * df/N,
        GCV =  logloss / (1 - df / N) ^ 2
      ))
    }
  }
