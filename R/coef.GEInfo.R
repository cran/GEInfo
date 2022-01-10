#' @title Extract coefficients from a fitted object
#' @description Report the estimate of all coefficients from a fitted "CGEInfo" or "GEInfo" model object.
#'
#' @param object A fitted "CGEInfo" or "GEInfo" model object for which the estimate of coefficients is extracted.
#' @param ... Other arguments.
#'
#' @return A coefficient vector of length (q+1) x (p+1), including the estimates for \eqn{\alpha} (intercept), \eqn{a} (coefficients for all E variables), and \eqn{b} (coefficients for all G variables and G-E interactions).

#' @export
coef.GEInfo <- function(object, ...)
{
  return(object$coef)
}
