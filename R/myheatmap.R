myheatmap <- function(G.coef, E.coef,  GE.coef){
  if(anyNA(names(G.coef)))
  {names(G.coef)<-paste0('G',1:length(G.coef))}
  if(anyNA(names(E.coef)))
  {names(E.coef)<-paste0('E',1:length(E.coef))}
  gname<-names(G.coef)
  ename<-names(E.coef)
  rownames(GE.coef) <- gname
  colnames(GE.coef) <- ename
  ann_row <- data.frame(G=G.coef)
  ann_col <- data.frame(E=E.coef)
  colors <- grDevices::heat.colors(100,rev=T)
  ann_colors <- list(G=colors,E=colors)
  pheatmap::pheatmap(GE.coef,annotation_row=ann_row,annotation_col=ann_col,
           annotation_colors=ann_colors,color=colors,
           cluster_cols=F, cluster_rows=F,
           fontsize_row=6,fontsize_col=6,fontsize = 6)
}
