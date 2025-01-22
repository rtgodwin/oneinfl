#' @title Expected Value for Zero-Truncated Negative Binomial Distribution
#' 
#' @description
#' Computes the expected value from a zero-truncated negative binomial (ZTNB) distribution 
#' based on the model parameters and covariates.
#' 
#' @param b Numeric vector of coefficients for the regression model.
#' @param a Numeric scalar, the overdispersion parameter of the negative binomial distribution.
#' @param X Matrix of predictors for the model, where rows correspond to observations 
#'   and columns to covariates.
#' 
#' @return
#' A numeric vector of expected values for each observation in the dataset.
#' 
#' @seealso
#' \code{\link{dEdq_nb_noinfl}} for computing partial derivatives of expected values in the ZTNB model.
#' 
#' @examples
#' \dontrun{
#' # Example inputs
#' b <- c(0.5, -0.2, 0.1)
#' a <- 1.5
#' X <- matrix(c(1, 2, 3, 1, 0.5, 1), ncol = 3)
#' 
#' # Compute expected values
#' E_negbin_noinfl(b, a, X)
#' }
#' 
#' @keywords internal

E_negbin_noinfl <- function(b, a, X) {
  l <- exp(X %*% b)
  th <- l / a
  w <- 0
  w + (1 - w) * (l / (1 - (1 + th) ^ (-a)))
}