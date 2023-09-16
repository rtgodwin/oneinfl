roipp <- function(b, g, X, Z) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  l <- exp(X %*% b)
  L <- -l / (exp(l) - l - 1)
  omega <- L + (1-L)/(1 + exp(-Z %*% g))
  n <- nrow(X)
  y <- rep(0, n)
  for(i in 1:n) {
    roll <- runif(1)
    # Determine prob of a 1 count
    probs <- omega[i] + (1 - omega[i]) * l / (exp(l) - 1)
    k <- 1
    # Determine subsequent probs until roll is reached
    while(probs[k] < roll) {
      k <- k + 1
      probs <- c(probs, probs[k - 1] + (1 - omega[i]) * l ^ k / ((exp(l) - 1) * factorial(k)))
    }
    y[i] <- k
  }
  return(y)
}
