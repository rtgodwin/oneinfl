#' @title Compute Partial Derivatives of Expected Values for a One-inflated Zero-truncated Negative Binomial Model
#' 
#' @description
#' This internal function computes the partial derivatives of expected values 
#' for a one-inflated zero-truncated negative binomial regression model with respect to covariates. 
#' It also accounts for marginal effects of dummy variables.
#' 
#' @param b Numeric vector of coefficients for the rate parameter.
#' @param g Numeric vector of coefficients for the inflation process.
#' @param a Numeric scalar, the overdispersion parameter of the negative binomial distribution.
#' @param X Matrix of predictors for the main model, where rows correspond to observations 
#'   and columns to covariates.
#' @param Z Matrix of predictors for the inflation process, structured similarly to `X`.
#' @param dummies Character vector of column names from `X` and `Z` that are considered 
#'   dummy variables for which marginal effects need to be computed.
#' 
#' @return
#' A matrix of partial derivatives (or marginal effects) with rows corresponding to observations 
#' and columns to covariates. For dummy variables, marginal effects are calculated directly.
#' 
#' @details
#' This function performs the following tasks:
#' - Computes partial derivatives of expected values with respect to covariates in `X` and `Z`.
#' - Handles marginal effects for dummy variables by comparing expected values when the dummy 
#'   variable is set to 0 versus 1.
#' 
#' The function is designed for internal use in the package and assumes that all input matrices 
#' and vectors are correctly specified. Any unexpected input structure may result in errors.
#' 
#' @seealso
#' \code{\link{E_negbin}} for expected values in the negative binomial model.
#' 
#' @keywords internal

dEdq_nb <- function(b, g, a, X, Z, dummies, formula) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  th <- l / a
  P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
  L <- -(((a / (a + l)) ^ (-a)) * (1 / l) * (1 + l/a - (1 + l/a) ^ (1-a)) - 1) ^ (-1)
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
  
  # Do the math
  dfdq <- dldq * as.vector(((1 + l / a) ^ (-a)) * (1 + l / a - (1 + l / a) ^ (1 - a)) ^ (-1) * (1 - (a ^ 2) * l / ((a + l) ^ 2) - (l / a) * (1 + l / a - (1 + l / a) ^ (1 - a)) ^ (-1) * (1 - (1 - a) * (1 + l / a) ^ (-a))))
  dLdq <- -dfdq / as.vector(((1 - P1) ^ 2))
  dwdq <- dLdq * as.vector((1 - 1 / (1 + t))) - dzdq * as.vector(((1 - P1) / (1 + t) ^ 2))
  dEdq <- dwdq * as.vector((1 - l / (1 - (1 + l / a) ^ (-a)))) + dldq * as.vector(((1 - w) / (1 - (1 + l / a) ^ (-a))) * (1 - (l * (1 + l / a) ^ (-a - 1)) / (1 - (1 + l / a) ^ (-a))))
  
  # Calculate the "marginal" effect for dummies properly
  if(length(dummies) == 0) {return(dEdq)}
  else {
    for(i in 1:length(dummies)) {
    Xd1 <- Xd0 <- X
    Zd1 <- Zd0 <- Z
    Xd1[ , dummies[i] == colnames(X)] <- 1
    Xd0[ , dummies[i] == colnames(X)] <- 0
    Zd1[ , dummies[i] == colnames(Z)] <- 1
    Zd0[ , dummies[i] == colnames(Z)] <- 0
    dEdq[, dummies[i]] <- E_negbin(b, g, a, X=Xd1, Z=Zd1) - E_negbin(b, g, a, X=Xd0, Z=Zd0)
    }
  return(dEdq)
  }
}