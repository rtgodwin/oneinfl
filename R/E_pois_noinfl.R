E_pois_noinfl <- function(b, X) {
  l <- exp(X %*% b)
  w <- 0
  w + (1 - w) * l * exp(l) / (exp(l) - 1)
}
