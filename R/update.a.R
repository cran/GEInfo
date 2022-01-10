update.a <- function(Eols, r) {
  fit <-stats::lm(r ~ Eols)
  if (length(unique(Eols[, 1])) == 1) {
    a0 <- stats::coef(fit)[-2]
  }  else {
    a0 <- stats::coef(fit)[-1]
  }
  return (list(a0 = a0))
}
