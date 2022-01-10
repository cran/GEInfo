#' @title Construct Type_S prior information
#' @description For G variables and G-E interactions, transform their prior information from counts(frequencies) into a set of significant variables (Type_S)
#' @param G.count A numeric vector, including the prior counts (frequencies) for G variables.
#' @param GE.count A numeric matrix, including the prior counts (frequencies) for G-E interactions.
#' @param eta_G A probability. The (eta_G)th quantile of G.count is used as a count (frequency) threshold (denoted by varphi_G) for G variables.
#'               Default is 0.95.
#' @param eta_GE A probability. The (eta_GE)th quantile of GE.count is used as a data-dependent count (frequency) threshold (denoted by varphi_GE) for G-E interactions.
#'              Default is 0.95.
#' @param varphi_G A user supplied  count threshold for G variables. It is used to determine which G variables will be finally included in the Type_S prior information set. Default is NULL.
#'               Typical usage is to have the program calculate the (eta_G)th quantile of G.count as the threshold. Supplying a varphi_G value will override this.
#' @param varphi_GE A user supplied threshold value used for G-E interactions. It is used to determine which G-E interactions will be finally included in the Type_S prior information set. Default is NULL.
#'                 Typical usage is to have the program calculate the (eta_GE)th quantile of GE.count as the threshold. Supplying a varphi_GE value will override this.
#' @return The outputs include the Type_S prior information sets for G variable and G-E interactions.
#'        \item{S_G}{A numeric vector, denoting the Type_S set for G variables. For j in S_G, the jth G variable
#'        is suggested to be associated with the response.}
#'       \item{S_GE}{A numeric matrix, denoting the Type_S set for G-E interactions. For (l,k) in S_GE,the lth G variable
#'        and the kth E variable is suggested to have an interaction effect on the response.}
#'
#' @export
#'
#' @examples
#' G.count<-c(100,300)
#' GE.count<-matrix(c(130,356,8,30,87,2),nrow=2)
#' TypeS(G.count,GE.count)
TypeS<-function(G.count,GE.count,eta_G=0.95, eta_GE=0.95, varphi_G=NULL, varphi_GE=NULL)
{
  if(!missing(GE.count))
  {
    if(is.null(varphi_GE)) varphi_GE<-stats::quantile(as.vector(GE.count),eta_GE)
    index.GE<-which(GE.count>=varphi_GE,arr.ind=TRUE)
    S_GE <-index.GE
  } else message('GE.count for GE interactions is missing')
  if(!missing(G.count))
  {
    if(is.null(varphi_G)) varphi_G<-stats::quantile(G.count,eta_G)
    index.G<-which(G.count>=varphi_G)
    S_G <-index.G
  }  else message('GE.count for GE interactions is missing')
  return(list(S_G=S_G, S_GE=S_GE))
}



