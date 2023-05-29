dfee_pois <- function(b, g, X, Z, dummies) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  L <- -l / (exp(l) - l - 1)
  w <- L + (1 - L) / (1 + t)
  
  dldq <- dzdq <- matrix(, nrow(X), (ncol(data) - 1))
  colnames(dldq) <- colnames(dzdq) <- colnames(data[-1])
  dldq[, colnames(X)[-1]] <- l %*% b[-1]
  dldq[, -which(colnames(dldq) %in% colnames(X)[-1])] <- 0L
  dzdq[, colnames(Z)[-1]] <- -t %*% g[-1]
  dzdq[, -which(colnames(dzdq) %in% colnames(Z)[-1])] <- 0L
  
  dwdq <- dldq * as.vector(((exp(l) - l * exp(l) - 1) / (exp(l) - l - 1) ^ 2) * ((-t) / (1 + t))) - dzdq * as.vector(((exp(l) - 1) / ((exp(l) - l - 1) * (1 + t) ^ 2)))
  dfeedq <- dwdq * as.vector((1 - (l * exp(l)) / (exp(l) - 1))) + dldq * as.vector(((w * exp(l) * (l - exp(l) + 1)) / ((exp(l) - 1) ^ 2)))
  
  for(i in 1:length(dummies)) {
    Xd1 <- Xd0 <- X
    Zd1 <- Zd0 <- Z
    Xd1[ , dummies[i] == colnames(X)] <- 1
    Xd0[ , dummies[i] == colnames(X)] <- 0
    Zd1[ , dummies[i] == colnames(Z)] <- 1
    Zd0[ , dummies[i] == colnames(Z)] <- 0
    dfeedq[, dummies[i]] <- fee_pois(b, g, X=Xd1, Z=Zd1) - fee_pois(b, g, X=Xd0, Z=Zd0)
  }
  dfeedq
}
