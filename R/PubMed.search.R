#' @title Search prior counts for G variables and G-E interactions
#' @description Provide an available tool for mining prior counts for G variables and G-E interactions from PubMed database.
#' @param Yname A user supplied character including disease name such as "breast".
#' @param Gname A user supplied character vector including all G variable names.
#' @param Ename A user supplied character vector including all E variable names.
#' @param Gnamefile A newline-delimited text file uploaded by users that contains all the G variable names to be searched. Each row represents a G variable name.
#'           It provides another way to input G variable names besides from argument "Gname".

#' @return Return the searched frequencies.
#' \item{G.count}{A numeric vector, presenting the prior counts for all searched G variables.}
#' \item{GE.count}{A numeric matrix of dimensions length(Gname) x length(Ename), which presents the prior counts
#'          for all G variables (Gname) and E variables (Ename) comparisons}
#' @export
#' @examples
#' Yname <- c('breast')
#' Gname <- c('CAMP')
#' Ename <- c('Age')
#' res <- PubMed.search(Yname,Gname,Ename)
#' res
#' @importFrom dplyr %>%
PubMed.search <- function(Yname,Gname,Ename, Gnamefile){
  requireNamespace('dplyr')
  if(!missing(Gnamefile)) Gname <- readLines(Gnamefile)
  G.count <- rep(0,length(Gname))
  names(G.count) <- Gname
  for(i in 1:length(Gname)){
    url <- paste0('https://pubmed.ncbi.nlm.nih.gov/?term=',Yname,'%20AND%20',Gname[i])
    html <- gethtml(url)
    num_res <- rvest::html_elements(html,'#search-results > div.top-wrapper > div.results-amount-container > div.results-amount > span')
    if(length(num_res)==0){
      G.count[i]<-0
    }else{
      count <- rvest::html_text2(num_res)
      count <- gsub(',','',count)
      count <- as.integer(count)
      G.count[i] <- count
    }
  }
  if(!missing(Ename)){
    GE.count <- matrix(nrow=length(Gname),ncol=length(Ename))
    rownames(GE.count) <- Gname
    colnames(GE.count) <- Ename
    for(i in 1:length(Gname)){
      for(j in 1:length(Ename)){
        url <- paste0('https://pubmed.ncbi.nlm.nih.gov/?term=',Yname,'%20AND%20',Gname[i],'%20AND%20',Ename[j])
        html <- gethtml(url)
        num_res <- rvest::html_elements(html,'#search-results > div.top-wrapper > div.results-amount-container > div.results-amount > span')
        if(length(num_res)==0){
          GE.count[i,j]<-0
        }else{
          count <- rvest::html_text2(num_res)
          count <- gsub(',','',count)
          count <- as.integer(count)
          GE.count[i,j] <- count
        }
      }
    }
  }
  if(missing(Ename)){
    res<-G.count
  }else{
    res<-list(G.count=G.count,GE.count=GE.count)
  }
  class(res) = 'PubMed'
  return(res)
}


