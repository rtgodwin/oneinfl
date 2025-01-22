#' @title Expected Value for One-Inflated Zero-Truncated Negative Binomial Distribution
#' 
#' @description
#' Computes the expected value from a one-inflated, zero-truncated negative binomial (OIZTNB) distribution. 
#' This function calculates the expected value based on the model parameters and covariates.
#' 
#' @param b Numeric vector of coefficients for the main model.
#' @param g Numeric vector of coefficients for the one-inflation process.
#' @param a Numeric scalar, the overdispersion parameter of the negative binomial distribution.
#' @param X Matrix of predictors for the main model, where rows correspond to observations 
#'   and columns to covariates.
#' @param Z Matrix of predictors for the one-inflation process, structured similarly to `X`.
#' 
#' @return
#' A numeric vector of expected values for each observation in the dataset.
#' 
#' @seealso
#' \code{\link{dEdq_nb}} for computing partial derivatives of expected values in the OIZTNB model.
#' 
#' @examples
#' \dontrun{
#' # Example inputs
#' b <- c(0.5, -0.2, 0.1)
#' g <- c(-1.2, 0.3)
#' a <- 1.5
#' X <- matrix(c(1, 2, 3, 1, 0.5, 1), ncol = 3)
#' Z <- matrix(c(1, 0.8, 1.5, 1, 1.2, 1), ncol = 2)
#' 
#' # Compute expected values
#' E_negbin(b, g, a, X, Z)
#' }
#' 
#' @keywords internal

E_negbin <- function(b, g, a, X, Z) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  th <- l / a
  P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
  L <- -P1 / (1 - P1)
  w <- L + (1 - L) / (1 + t)
  w + (1 - w) * (l / (1 - (1 + th) ^ (-a)))
}