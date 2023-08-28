roiztnb <- function(b, g, alpha, X, Z) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  l <- exp(X %*% b)
  L <- -(((alpha / (alpha + l)) ^ (-alpha)) * (1 / l) * (1 + l/alpha - (1 + l/alpha) ^ (1-alpha)) - 1) ^ (-1)
  omega <- L + (1-L)/(1 + exp(-Z %*% g))
  theta <- l / alpha
  n <- nrow(X)
  y <- rep(0, n)
  for(i in 1:n) {
    roll <- runif(1)
    probs <- omega[i] + (1 - omega[i]) * alpha * ((1 / (1 + theta[i])) ^ alpha) * (theta[i] / (1 + theta[i] - (1 + theta[i]) ^ (1 - alpha)))
    k <- 1
    while(probs[k] < roll) {
      k <- k + 1
      probs <- c(probs, probs[k - 1] + (1 - omega[i]) * ((gamma(alpha + k)) / (gamma(alpha) * gamma(k + 1))) * ((1/(1 + theta[i])) ^ alpha) * ((theta[i] / (1 + theta[i])) ^ k) * (1 / (1 - (1 + theta[i]) ^ (-alpha))))
    }
    y[i] <- k
  }
  return(y)
}
