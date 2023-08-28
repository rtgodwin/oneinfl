dEdq_nb <- function(b, g, a, X, Z, dummies) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  th <- l / a
  P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
  L <- -(((a / (a + li)) ^ (-a)) * (1 / l) * (1 + l/a - (1 + l/a) ^ (1-a)) - 1) ^ (-1)
  w <- L + (1 - L) / (1 + t)
  
  dldq <- dzdq <- matrix(, nrow(X), (ncol(data) - 1))
  colnames(dldq) <- colnames(dzdq) <- colnames(data[-1])
  dldq[, colnames(X)[-1]] <- l %*% b[-1]
  dldq[, -which(colnames(dldq) %in% colnames(X)[-1])] <- 0L
  dzdq[, colnames(Z)[-1]] <- -t %*% g[-1]
  dzdq[, -which(colnames(dzdq) %in% colnames(Z)[-1])] <- 0L
  
  # Do the math
  dfdq <- dldq * as.vector(((1 + l / a) ^ (-a)) * (1 + l / a - (1 + l / a) ^ (1 - a)) ^ (-1) * (1 - (a ^ 2) * l / ((a + l) ^ 2) - (l / a) * (1 + l / a - (1 + l / a) ^ (1 - a)) ^ (-1) * (1 - (1 - a) * (1 + l / a) ^ (-a))))
  dLdq <- -dfdq / as.vector(((1 - P1) ^ 2))
  dwdq <- dLdq * as.vector((1 - 1 / (1 + t))) - dzdq * as.vector(((1 - P1) / (1 + t) ^ 2))
  dEdq <- dwdq * as.vector((1 - l / (1 - (1 + l / a) ^ (-a)))) + dldq * as.vector(((1 - w) / (1 - (1 + l / a) ^ (-a))) * (1 - (l * (1 + l / a) ^ (-a - 1)) / (1 - (1 + l / a) ^ (-a))))
  
  # Calculate the "marginal" effect for dummies properly
  if(length(dummies) == 0) {return(dEdq)}
  else {
    for(i in 1:length(dummies)) {
    Xd1 <- Xd0 <- X
    Zd1 <- Zd0 <- Z
    Xd1[ , dummies[i] == colnames(X)] <- 1
    Xd0[ , dummies[i] == colnames(X)] <- 0
    Zd1[ , dummies[i] == colnames(Z)] <- 1
    Zd0[ , dummies[i] == colnames(Z)] <- 0
    dEdq[, dummies[i]] <- E_negbin(b, g, a, X=Xd1, Z=Zd1) - E_negbin(b, g, a, X=Xd0, Z=Zd0)
    }
  return(dEdq)
  }
}
