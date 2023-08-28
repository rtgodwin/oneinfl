E_pois <- function(b, g, X, Z) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  w <- -l * (exp(l) - l - 1) ^ -1 + (1 + l * (exp(l) - l - 1) ^ -1) * (1 + t) ^ -1
  w * (1 - l * exp(l) / (exp(l) - 1))
}
