#' Generate Random Counts from a One-Inflated Poisson Process
#'
#' Simulates count data from a one-inflated Poisson process using specified parameters 
#' for the rate and one-inflation components.
#'
#' @param b A numeric vector of coefficients for the rate component.
#' @param g A numeric vector of coefficients for the one-inflation component.
#' @param X A matrix or data frame of predictor variables for the rate component.
#' @param Z A matrix or data frame of predictor variables for the one-inflation component.
#'
#' @return
#' A numeric vector of simulated count data.
#'
#' @details
#' This function generates count data from a one-inflated Poisson process. The process 
#' combines:
#' \itemize{
#'   \item A Poisson distribution for counts greater than one.
#'   \item A one-inflation component that adjusts the probability of observing a count of one.
#' }
#'
#' The algorithm:
#' \enumerate{
#'   \item Calculates the rate parameter (\eqn{\lambda}) as \eqn{\exp(X \cdot b)}.
#'   \item Computes the one-inflation probabilities (\eqn{\omega}) based on \eqn{Z \cdot g}.
#'   \item Simulates counts for each observation:
#'     \itemize{
#'       \item Draws a random number to determine whether the count is one.
#'       \item Iteratively calculates probabilities for higher counts until the random number is matched.
#'     }
#' }
#'
#' This function is useful for generating synthetic data for testing or simulation studies involving 
#' one-inflated Poisson models.
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
#' simulated_data <- roipp(b, g, X, Z)
#' print(simulated_data)
#'
#' @export

roipp <- function(b, g, X, Z) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  l <- exp(X %*% b)
  L <- -l / (exp(l) - l - 1)
  omega <- L + (1 - L) / (1 + exp(-Z %*% g))
  n <- nrow(X)
  y <- rep(0, n)
  for (i in 1:n) {
    roll <- runif(1)
    # Determine prob of a 1 count
    probs <- omega[i] + (1 - omega[i]) * l[i] / (exp(l[i]) - 1)
    k <- 1
    # Determine subsequent probs until roll is reached
    while (probs[k] < roll) {
      k <- k + 1
      probs <- c(probs, probs[k - 1] + (1 - omega[i]) * (l[i] ^ k) / ((exp(l[i]) - 1) * factorial(k)))
    }
    y[i] <- k
  }
  return(y)
}