#' @title Expected Value for Positive Poisson Distribution
#' 
#' @description
#' Computes the expected value from a Poisson (PP) distribution based on the model parameters and covariates.
#' 
#' @param b Numeric vector of coefficients for the Poisson model.
#' @param X Matrix of predictors for the Poisson model, where rows correspond 
#'   to observations and columns to covariates.
#' 
#' @return
#' A numeric vector of expected values for each observation in the dataset.
#' 
#' @seealso
#' \code{\link{dEdq_pois_noinfl}} for computing partial derivatives of expected values in the Poisson model.
#' 
#' @keywords internal

E_pois_noinfl <- function(b, X) {
  l <- exp(X %*% b)
  w <- 0
  w + (1 - w) * l * exp(l) / (exp(l) - 1)
}