#' @title Compute Partial Derivatives of Expected Values for Zero-truncated Negative Binomial Model
#' 
#' @description
#' This internal function calculates the partial derivatives of expected values 
#' for a regular truncated negative binomial regression model with respect to covariates. 
#' It also adjusts for marginal effects of dummy variables when specified.
#' 
#' @param b Numeric vector of coefficients for the regression model.
#' @param a Numeric scalar, the overdispersion parameter of the negative binomial distribution.
#' @param X Matrix of predictors for the model, where rows correspond to observations 
#'   and columns to covariates.
#' @param dummies Character vector of column names from `X` that are considered 
#'   dummy variables for which marginal effects need to be computed.
#' 
#' @return
#' A matrix of partial derivatives (or marginal effects) with rows corresponding to observations 
#' and columns to covariates. Marginal effects for dummy variables are calculated by contrasting 
#' expected values when the dummy variable is set to 0 versus 1.
#' 
#' @details
#' This function performs the following tasks:
#' - Computes partial derivatives of expected values with respect to covariates in `X`.
#' - Adjusts for marginal effects of dummy variables by altering their values in the design matrix 
#'   and computing the difference in expected values.
#' 
#' It is designed for internal use and assumes correct input structure. Improper inputs may result 
#' in errors or unexpected behavior.
#' 
#' @seealso
#' \code{\link{E_negbin_noinfl}} for computing expected values in the regular truncated negative binomial model.
#' 
#' @keywords internal
#' 

dEdq_nb_noinfl <- function(b, a, X, dummies, formula) {
  l <- exp(X %*% b)
  th <- l / a
  w <- 0
  
  vars <- attr(terms(formula), "term.labels")
  num_unique <- length(unique(vars))
  
  dldq <- matrix(, nrow(X), num_unique)
  colnames(dldq) <- vars
  dldq[, colnames(X)[-1]] <- l %*% b[-1]
  dldq[, -which(colnames(dldq) %in% colnames(X)[-1])] <- 0L

  # Do the math
  dwdq <- 0
  dEdq <- dwdq * as.vector((1 - l / (1 - (1 + l / a) ^ (-a)))) + dldq * as.vector(((1 - w) / (1 - (1 + l / a) ^ (-a))) * (1 - (l * (1 + l / a) ^ (-a - 1)) / (1 - (1 + l / a) ^ (-a))))
  
  # Calculate the "marginal" effect for dummies properly
  if(length(dummies) == 0) {return(dEdq)}
  else {
    for(i in 1:length(dummies)) {
      Xd1 <- Xd0 <- X
      Xd1[ , dummies[i] == colnames(X)] <- 1
      Xd0[ , dummies[i] == colnames(X)] <- 0
      dEdq[, dummies[i]] <- E_negbin_noinfl(b, a, X=Xd1) - E_negbin_noinfl(b, a, X=Xd0)
    }
    return(dEdq)
  }
}