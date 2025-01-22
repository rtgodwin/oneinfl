#' @title Expected Value for One-Inflated Poisson Distribution
#' 
#' @description
#' Computes the expected value from a one-inflated Poisson (OIPP) distribution 
#' based on the model parameters and covariates.
#' 
#' @param b Numeric vector of coefficients for the main Poisson model.
#' @param g Numeric vector of coefficients for the one-inflation process.
#' @param X Matrix of predictors for the main Poisson model, where rows correspond 
#'   to observations and columns to covariates.
#' @param Z Matrix of predictors for the one-inflation process, structured similarly to `X`.
#' 
#' @return
#' A numeric vector of expected values for each observation in the dataset.
#' 
#' @seealso
#' \code{\link{dEdq_pois}} for computing partial derivatives of expected values in the OIPP model.
#' 
#' @examples
#' \dontrun{
#' # Example inputs
#' b <- c(0.5, -0.2, 0.1)
#' g <- c(-0.3, 0.7)
#' X <- matrix(c(1, 2, 3, 1, 0.5, 1), ncol = 3)
#' Z <- matrix(c(1, 0.8, 1.5, 1, 1.2, 1), ncol = 2)
#' 
#' # Compute expected values
#' E_pois(b, g, X, Z)
#' }
#' 
#' @keywords internal

E_pois <- function(b, g, X, Z) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  w <- -l * (exp(l) - l - 1) ^ -1 + (1 + l * (exp(l) - l - 1) ^ -1) * (1 + t) ^ -1
  w + (1 - w) * l * exp(l) / (exp(l) - 1)
}