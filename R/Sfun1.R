Sfun1 <- function(z, c) {
  value <- max((1 - c / abs(z)), 0) * z
  return(value)
}
