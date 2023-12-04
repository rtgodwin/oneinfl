# Caluculate the expected value from an OIZTNB distribution
E_negbin <- function(b, g, a, X, Z) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  th <- l / a
  P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
  L <- -P1 / (1 - P1)
  w <- L + (1 - L) / (1 + t)
  w + (1 - w) * (l / (1 - (1 + th) ^ (-a)))
}
