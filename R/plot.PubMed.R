#' @title   Visualize the prior information
#' @description  Visualize the prior counts for G variables and G-E interactions. It reports a bar chart for the top 40 G variables by prior count and
#'        a boxplot of prior counts for all G variables. For each E variables, it draws a bar chart for the corresponding top 20 G-E interactions by prior count.

#' @param x A 'PubMed' object for which visualization is desired.
#' @param ... Other parameters.
#' @param G.count A numeric vector of length p, including prior counts for all G variables. Default is NULL.
#' @param GE.count A numeric matrix of dimensions p*q, including prior counts for G-E interactions. Default is NULL.

#' @return The output includes bar chart for top G variables and G-E interactions by prior counts, and boxplot of prior counts for all G variables.
#' @export
#' @importFrom graphics title text layout barplot boxplot
#' @importFrom grDevices rainbow
plot.PubMed <-function (x,G.count=NULL,GE.count=NULL,...)
{
if(is.null(G.count))
 { G.count<-x$G.count}
if(is.null(GE.count))
{GE.count<-x$GE.count}

nonzero_g <- G.count[which(G.count!=0)]
max_idx <- ifelse(length(nonzero_g)<=20,length(nonzero_g),20)
top20_g <- sort(nonzero_g,decreasing=T)[1:max_idx]
graphics::layout(matrix(c(1,2),1,2),widths=c(2,1))
pic<-graphics::barplot(top20_g,col=grDevices::rainbow(length(top20_g)),space=1.2,
                       axisnames=T,names.arg=names(top20_g),
                       axis.lty = 0.4,las=2,cex.names = 0.45,
                       cex.lab=0.7,cex.axis=0.5,horiz=T,xlab='Prior count',
                       xlim=c(0,top20_g[1]*1.2+1))
graphics::text(y=pic,x=top20_g,labels=top20_g,cex=0.4,pos=4)
graphics:: boxplot(G.count,cex.axis=0.7,color='red',horizontal = F,
                   ylab='Prior count',cex.lab=0.7,pars=list(outwex=0.2))

if(!missing(GE.count)){
  for(i in 1:ncol(GE.count)){
    if(i%%4==1) graphics::layout(matrix(c(1,2,3,4),2,2,byrow=T))
    nonzero_ge <- GE.count[,i][which(GE.count[,i]!=0)]
    max_idx <- ifelse(length(nonzero_ge)<=20,length(nonzero_ge),20)
    top20_ge <- sort(nonzero_ge,decreasing=T)[1:max_idx]
    pic<-graphics::barplot(top20_ge,col=rainbow(length(top20_ge)),space=1.2,
                           axisnames=T,names.arg=names(top20_ge),
                           axis.lty = 0.4,las=2,cex.names = 0.45,
                           cex.lab=0.7,cex.axis=0.5,horiz=T,xlab=colnames(GE.count) [i],
                           xlim=c(0,top20_ge[1]*1.2+1))
    graphics::text(y=pic,x=top20_ge,labels=top20_ge,cex=0.4,pos=4)
  }
}
}
