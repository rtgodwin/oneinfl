rpp <- function(b, X) {
  X <- as.matrix(X)
  l <- exp(X %*% b)
  n <- nrow(X)
  y <- rep(0, n)
  for(i in 1:n) {
    roll <- runif(1)
    # Determine prob of a 1 count
    probs <- l[i] / (exp(l[i]) - 1)
    k <- 1
    # Determine subsequent probs until roll is reached
    while(probs[k] < roll) {
      k <- k + 1
      probs <- c(probs, probs[k - 1] + (l[i] ^ k) / ((exp(l[i]) - 1) * factorial(k)))
    }
    y[i] <- k
  }
  return(y)
}
