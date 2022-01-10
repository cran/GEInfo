#' @title Make Predictions for a fitted model
#' @description  Output predicted response values for new observations.
#' @param Enew Matrix of dimensions \eqn{n_test x q} for E variables at which predictions are to be made.
#' @param Gnew Matrix of dimensions \eqn{n_test x p}for G variables at which predictions are to be made.
#' @param family Model type: one of ("gaussian", "binomial", "poisson").
#' @param object A fitted "GEInfo" model object for which prediction is desired.
#' @param ... Other arguments.
#'
#' @return Return a vector of length \eqn{n_test}, representing the fitted response value. For family= “gaussian”, the fitted values are returned;
#' for family = “binary”, the fitted probabilities are returned;
#' for family = “poisson”, the fitted means are returned.
#' @export
predict.GEInfo <- function(object, Enew, Gnew, family, ...)
{
  Wnew <- matW(Enew, Gnew)
  coef.res <- object$coef
  eta.new <- cbind(1, Enew, Wnew) %*% coef.res
  pred.Y <- switch (
    family,
    gaussian = eta.new,
    binomial = 1 / (1 + exp (-eta.new)),
    poisson = exp (eta.new)
  )
  return (pred.Y)
}
