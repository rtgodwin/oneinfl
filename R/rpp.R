#' Generate Random Counts from a Zero-Truncated Poisson Process
#'
#' Simulates count data from a zero-truncated Poisson process using specified parameters 
#' for the rate component.
#'
#' @param b A numeric vector of coefficients for the rate component.
#' @param X A matrix or data frame of predictor variables for the rate component.
#'
#' @return
#' A numeric vector of simulated count data.
#'
#' @details
#' This function generates count data from a zero-truncated Poisson process, which models 
#' count data without zeros. The process involves:
#' \itemize{
#'   \item Calculating the rate parameter (\eqn{\lambda}) as \eqn{\exp(X \cdot b)}.
#'   \item Iteratively computing probabilities for counts starting from 1 and adding to the cumulative 
#'     probability until a randomly drawn value is matched.
#' }
#'
#' This function is useful for generating synthetic data for testing or simulation studies involving 
#' zero-truncated Poisson models.
#'
#' @seealso
#' \code{\link{truncreg}} for fitting zero-truncated models.
#'
#' @examples
#' # Example usage
#' set.seed(123)
#' X <- matrix(rnorm(100), ncol = 2)
#' b <- c(0.5, -0.2)
#' simulated_data <- rpp(b, X)
#' print(simulated_data)
#'
#' @export
rpp <- function(b, X) {
  X <- as.matrix(X)
  l <- exp(X %*% b)
  n <- nrow(X)
  y <- rep(0, n)
  for (i in 1:n) {
    roll <- runif(1)
    # Determine prob of a 1 count
    probs <- l[i] / (exp(l[i]) - 1)
    k <- 1
    # Determine subsequent probs until roll is reached
    while (probs[k] < roll) {
      k <- k + 1
      probs <- c(probs, probs[k - 1] + (l[i] ^ k) / ((exp(l[i]) - 1) * factorial(k)))
    }
    y[i] <- k
  }
  return(y)
}