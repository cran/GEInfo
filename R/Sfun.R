Sfun <- function(z,c) {
  znorm <- sqrt(sum(z^2))
  value <- max((1-c/znorm),0)*z
  return (value)
}
