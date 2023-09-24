dEdq_nb_noinfl <- function(b, a, X, dummies) {
  l <- exp(X %*% b)
  th <- l / a
  w <- 0
  
  dldq <- matrix(, nrow(X), (ncol(data) - 1))
  colnames(dldq) <- colnames(data[-1])
  dldq[, colnames(X)[-1]] <- l %*% b[-1]
  dldq[, -which(colnames(dldq) %in% colnames(X)[-1])] <- 0L

  # Do the math
  dwdq <- 0
  dEdq <- dwdq * as.vector((1 - l / (1 - (1 + l / a) ^ (-a)))) + dldq * as.vector(((1 - w) / (1 - (1 + l / a) ^ (-a))) * (1 - (l * (1 + l / a) ^ (-a - 1)) / (1 - (1 + l / a) ^ (-a))))
  
  # Calculate the "marginal" effect for dummies properly
  if(length(dummies) == 0) {return(dEdq)}
  else {
    for(i in 1:length(dummies)) {
      Xd1 <- Xd0 <- X
      Xd1[ , dummies[i] == colnames(X)] <- 1
      Xd0[ , dummies[i] == colnames(X)] <- 0
      dEdq[, dummies[i]] <- E_negbin_noinfl(b, a, X=Xd1) - E_negbin_noinfl(b, a, X=Xd0)
    }
    return(dEdq)
  }
}
