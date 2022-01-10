gethtml <- function(url){
  repeat{
    ret <- try(rvest::read_html(url),silent=T)
    if('try-error' %in% class(ret)){
      next
    }else{
      html <- ret
      break
    }
  }
  return(html)
}
