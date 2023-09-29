dEdq_pois_noinfl <- function(b, X, dummies) {
  l <- exp(X %*% b)
  w <- 0
  
  dldq <- matrix(, nrow(X), (ncol(data) - 1))
  colnames(dldq) <- colnames(data[-1])
  dldq[, colnames(X)[-1]] <- l %*% b[-1]
  dldq[, -which(colnames(dldq) %in% colnames(X)[-1])] <- 0L

  # Do the math
  dwdq <- 0
  dEdq <- dwdq * as.vector((1 - (l * exp(l)) / (exp(l) - 1))) + dldq * as.vector((1 - w) * exp(l) * (exp(l) - l - 1) / (exp(l) - 1) ^ 2)
  
  # Calculate the "marginal" effect for dummies properly
  if(length(dummies) == 0) {return(dEdq)}
  else {
    for(i in 1:length(dummies)) {
      Xd1 <- Xd0 <- X
      Xd1[ , dummies[i] == colnames(X)] <- 1
      Xd0[ , dummies[i] == colnames(X)] <- 0
      dEdq[, dummies[i]] <- E_pois_noinfl(b, X=Xd1) - E_pois_noinfl(b, X=Xd0)
    }
    return(dEdq)
  }
}
