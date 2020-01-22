dfee_nb <- function(b, g, a, X, Z, dummies) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  th <- l / a
  P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
  L <- -P1 / (1 - P1)
  w <- L + (1 - L) / (1 + t)

  dldq <- dzdq <- matrix(, nrow(X), (ncol(data) - 1))

  colnames(dldq) <- colnames(dzdq) <- colnames(data[-1])

  dldq[, colnames(X)[-1]] <- l %*% b[-1]
  dldq[, -which(colnames(dldq) %in% colnames(X)[-1])] <- 0L

  dzdq[, colnames(Z)[-1]] <- -t %*% g[-1]
  dzdq[, -which(colnames(dzdq) %in% colnames(Z)[-1])] <- 0L

  dfdq <- dldq * as.vector(((1 + l / a) ^ (-a)) * (1 + l / a - (1 + l / a) ^ (1 - a)) ^ (-1) * (1 - (a ^ 2) * l / ((a + l) ^ 2) - (l / a) * (1 + l / a - (1 + l / a) ^ (1 - a)) ^ (-1) * (1 - (1 - a) * (1 + l / a) ^ (-a))))

  dLdq <- -dfdq / as.vector(((1 - P1) ^ 2))

  dwdq <- dLdq * as.vector((1 - 1 / (1 + t))) - dzdq * as.vector(((1 - P1) / (1 + t) ^ 2))

  dfeedq <- dwdq * as.vector((1 - l / (1 - (1 + l / a) ^ (-a)))) - dldq * as.vector((w / (1 - (1 + l / a) ^ (-a))) * (1 + (l * (1 + l / a) ^ (-a - 1)) / (1 - (1 + l / a) ^ (-a))))

  dfeedq
}
