#' Generate Random Counts from a One-Inflated Zero-Truncated Negative Binomial Process
#'
#' Simulates count data from a one-inflated, zero-truncated negative binomial (OIZTNB) process 
#' using specified parameters for the rate, one-inflation, and dispersion components.
#'
#' @param b A numeric vector of coefficients for the rate component.
#' @param g A numeric vector of coefficients for the one-inflation component.
#' @param alpha A numeric value representing the dispersion parameter for the negative binomial distribution.
#' @param X A matrix or data frame of predictor variables for the rate component.
#' @param Z A matrix or data frame of predictor variables for the one-inflation component.
#'
#' @return
#' A numeric vector of simulated count data.
#'
#' @details
#' This function generates count data from a one-inflated, zero-truncated negative binomial process. 
#' The process combines:
#' \itemize{
#'   \item A negative binomial distribution for counts greater than one.
#'   \item A one-inflation component that adjusts the probability of observing a count of one.
#' }
#'
#' The algorithm:
#' \enumerate{
#'   \item Calculates the rate parameter (\eqn{\lambda}) as \eqn{\exp(X \cdot b)}.
#'   \item Computes the one-inflation probabilities (\eqn{\omega}) based on \eqn{Z \cdot g}.
#'   \item Computes the negative binomial dispersion parameter (\eqn{\theta = \lambda / \alpha}).
#'   \item Simulates counts for each observation:
#'     \itemize{
#'       \item Draws a random number to determine whether the count is one.
#'       \item Iteratively calculates probabilities for higher counts until the random number is matched.
#'     }
#' }
#'
#' This function is useful for generating synthetic data for testing or simulation studies involving 
#' one-inflated, zero-truncated negative binomial models.
#'
#' @seealso
#' \code{\link{oneinfl}} for fitting one-inflated models.
#'
#' @examples
#' # Example usage
#' set.seed(123)
#' X <- matrix(rnorm(100), ncol = 2)
#' Z <- matrix(rnorm(100), ncol = 2)
#' b <- c(0.5, -0.2)
#' g <- c(1.0, 0.3)
#' alpha <- 1.5
#' simulated_data <- roiztnb(b, g, alpha, X, Z)
#' print(simulated_data)
#'
#' @export

roiztnb <- function(b, g, alpha, X, Z) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  l <- exp(X %*% b)
  L <- -(((alpha / (alpha + l)) ^ (-alpha)) * (1 / l) * (1 + l / alpha - (1 + l / alpha) ^ (1 - alpha)) - 1) ^ (-1)
  omega <- L + (1 - L) / (1 + exp(-Z %*% g))
  theta <- l / alpha
  n <- nrow(X)
  y <- rep(0, n)
  for (i in 1:n) {
    roll <- runif(1)
    probs <- omega[i] + (1 - omega[i]) * alpha * ((1 / (1 + theta[i])) ^ alpha) * (theta[i] / (1 + theta[i] - (1 + theta[i]) ^ (1 - alpha)))
    k <- 1
    while (probs[k] < roll) {
      k <- k + 1
      probs <- c(probs, probs[k - 1] + (1 - omega[i]) * ((gamma(alpha + k)) / (gamma(alpha) * gamma(k + 1))) * ((1 / (1 + theta[i])) ^ alpha) * ((theta[i] / (1 + theta[i])) ^ k) * (1 / (1 - (1 + theta[i]) ^ (-alpha))))
      if (is.nan(probs[k]) == TRUE) { break }
    }
    y[i] <- k
  }
  return(y)
}
