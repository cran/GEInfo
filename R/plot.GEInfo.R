#' @title Heatmap of the identification results
#' @description Plot the heatmap for all E variables, identified G variables, and their G-E interactions from a fitted (GEInfo) model.
#'
#' @param x A fitted "GEInfo" model object for which prediction is desired.
#' @param Gname Names of all G variables. Default is NULL.
#' @param Ename Names of all E variables. Default is NULL.
#' @param ... Other parameters.
#' @return A Heatmap.

#' @export

plot.GEInfo<-function(x,Gname=NULL,Ename=NULL,...)
{
  G.coef0 <- x$beta
  selct.index <- which(G.coef0!=0)
  E.coef <- x$a
  if(!is.null(Gname)) names(G.coef0) <- Gname
  if(!is.null(Ename)) names(E.coef) <- Ename

  GE.coef0 <- x $gamma
  G.coef <- G.coef0[selct.index]
  GE.coef <-matrix(GE.coef0[selct.index,],ncol=length(E.coef),byrow=TRUE)
  myheatmap(G.coef, E.coef, GE.coef)
}
