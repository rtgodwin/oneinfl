dEdq_pois_noinfl <- function(b, X, dummies) {
  l <- exp(X %*% b)
    
  dldq <- matrix(, nrow(X), (ncol(X) - 1))
  colnames(dldq) <- colnames(X[-1])
  dldq[, colnames(X)[-1]] <- l %*% b[-1]
  dldq[, -which(colnames(dldq) %in% colnames(X)[-1])] <- 0L

  # Do the math
  dEdq <- dldq * exp(l) * (exp(l) - l - 1) / ((exp(l) - 1) ^ 2)
  
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
