#' @title Compute Partial Derivatives of Expected Values for a Positive Poisson Model
#' 
#' @description
#' This internal function calculates the partial derivatives of expected values 
#' for a positive Poisson regression model with respect to covariates. 
#' It also computes marginal effects for specified dummy variables.
#' 
#' @param b Numeric vector of coefficients for the Poisson model.
#' @param X Matrix of predictors for the Poisson model, where rows correspond 
#'   to observations and columns to covariates.
#' @param dummies Character vector of column names from `X` that are treated 
#'   as dummy variables for which marginal effects are computed.
#' 
#' @return
#' A matrix of partial derivatives (or marginal effects) with rows corresponding to observations 
#' and columns to covariates. Marginal effects for dummy variables are calculated by contrasting 
#' expected values when the dummy is set to 0 versus 1.
#' 
#' @details
#' This function:
#' - Computes partial derivatives of expected values with respect to covariates in `X`.
#' - Handles marginal effects for dummy variables by modifying their values in the design matrix 
#'   and computing the difference in expected values.
#' 
#' It is designed for internal use and assumes correct input structure. Improper inputs may result 
#' in errors or unexpected behavior.
#' 
#' @seealso
#' \code{\link{E_pois_noinfl}} for computing expected values in the positive Poisson model.
#' 
#' @keywords internal

dEdq_pois_noinfl <- function(b, X, dummies, formula) {
  l <- exp(X %*% b)
  
  vars <- attr(terms(formula), "term.labels")
  num_unique <- length(unique(vars))
    
  dldq <- matrix(, nrow(X), num_unique)
  colnames(dldq) <- vars
  dldq[, colnames(X)[-1]] <- l %*% b[-1]
  dldq[, -which(colnames(dldq) %in% colnames(X)[-1])] <- 0L

  # Do the math
  dEdq <- dldq * exp(l) * (exp(l) - l - 1) / ((exp(l) - 1) ^ 2)
  
  # Calculate the "marginal" effect for dummies properly
  if(length(dummies) == 0) {return(dEdq)}
  else {
    for(i in 1:length(dummies)) {
      Xd1 <- Xd0 <- X
      Xd1[ , dummies[i] == colnames(X)] <- 1
      Xd0[ , dummies[i] == colnames(X)] <- 0
      dEdq[, dummies[i]] <- E_pois_noinfl(b, X=Xd1) - E_pois_noinfl(b, X=Xd0)
    }
    return(dEdq)
  }
}