#' @title Compute Partial Derivatives of Expected Values for a One-Inflated Positive Poisson Model
#' 
#' @description
#' This internal function calculates the partial derivatives of expected values 
#' for a one-inflated Poisson regression model with respect to covariates. 
#' It also computes marginal effects for specified dummy variables.
#' 
#' @param b Numeric vector of coefficients for the main Poisson model.
#' @param g Numeric vector of coefficients for the one-inflation process.
#' @param X Matrix of predictors for the main Poisson model, where rows correspond 
#'   to observations and columns to covariates.
#' @param Z Matrix of predictors for the one-inflation process, structured similarly to `X`.
#' @param dummies Character vector of column names from `X` and `Z` that are treated 
#'   as dummy variables for which marginal effects are computed.
#' 
#' @return
#' A matrix of partial derivatives (or marginal effects) with rows corresponding to observations 
#' and columns to covariates. For dummy variables, marginal effects are calculated by contrasting 
#' expected values when the dummy is set to 0 versus 1.
#' 
#' @details
#' This function:
#' - Computes partial derivatives of expected values with respect to covariates in `X` and `Z`.
#' - Handles marginal effects for dummy variables by modifying their values in the design matrices 
#'   and computing the difference in expected values.
#' 
#' It is designed for internal use and assumes correct input structure. Improper inputs may result 
#' in errors or unexpected behavior.
#' 
#' @seealso
#' \code{\link{E_pois}} for computing expected values in the one-inflated Poisson model.
#' 
#' @keywords internal

dEdq_pois <- function(b, g, X, Z, dummies, formula) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  L <- -l / (exp(l) - l - 1)
  w <- L + (1 - L) / (1 + t)
  
  vars <- attr(terms(formula), "term.labels")
  expanded_vars <- unlist(strsplit(vars, "\\|"))
  expanded_vars <- trimws(expanded_vars)
  num_unique <- length(unique(expanded_vars))
  
  dldq <- dzdq <- matrix(, nrow(X), num_unique)
  colnames(dldq) <- colnames(dzdq) <- expanded_vars
  dldq[, colnames(X)[-1]] <- l %*% b[-1]
  dldq[, -which(colnames(dldq) %in% colnames(X)[-1])] <- 0L
  dzdq[, colnames(Z)[-1]] <- -t %*% g[-1]
  dzdq[, -which(colnames(dzdq) %in% colnames(Z)[-1])] <- 0L
  
  dwdq <- dldq * as.vector(((exp(l) - l * exp(l) - 1) / (exp(l) - l - 1) ^ 2) * ((-t) / (1 + t))) - dzdq * as.vector(((exp(l) - 1) / ((exp(l) - l - 1) * (1 + t) ^ 2)))
  dEdq <- dwdq * as.vector((1 - (l * exp(l)) / (exp(l) - 1))) + dldq * as.vector((1 - w) * exp(l) * (exp(l) - l - 1) / (exp(l) - 1) ^ 2)

  if(length(dummies) == 0) {return(dEdq)}
  else {
    for(i in 1:length(dummies)) {
      Xd1 <- Xd0 <- X
      Zd1 <- Zd0 <- Z
      Xd1[ , dummies[i] == colnames(X)] <- 1
      Xd0[ , dummies[i] == colnames(X)] <- 0
      Zd1[ , dummies[i] == colnames(Z)] <- 1
      Zd0[ , dummies[i] == colnames(Z)] <- 0
      dEdq[, dummies[i]] <- E_pois(b, g, X=Xd1, Z=Zd1) - E_pois(b, g, X=Xd0, Z=Zd0)
    }
    return(dEdq)
  }
}