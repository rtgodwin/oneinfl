E_negbin_noinfl <- function(b, a, X) {
  l <- exp(X %*% b)
  th <- l / a
  w <- 0
  w + (1 - w) * (l / (1 - (1 + th) ^ (-a)))
}
