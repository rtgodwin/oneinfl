#' @title Compute Marginal Effects for One-inflated models
#' 
#' @description
#' This wrapper function calls a different function to calculate marginal effects depending on the model type.
#' The marginal effects of the variables are evaluated at specified points, such as the sample means or averages, or at custom-defined cases.
#'
#' @param model An object representing a fitted model. Must be of class `oneinflmodel` or `truncmodel`.
#' @param df A \code{data.frame} containing the dataset used to fit the model. The variables in the data must match those used in the model.
#' @param at A character string or list specifying where to evaluate the marginal effects:
#' \itemize{
#'   \item \code{"AE"}: Average Effect (marginal effect averaged over all data points; default).
#'   \item \code{"EM"}: Effect at Means (marginal effect evaluated at the sample means of the data).
#'   \item A named list: A custom case specifying representative values for variables.
#' }
#' @param verbose Logical; if \code{TRUE} (default), prints the summary output. If \code{FALSE}, suppresses output table and returns a list containing several components.
#'
#' @return If \code{verbose=TRUE} (default), prints the marginal effects, their standard errors, z-values, p-values, and significance levels.
#' @return If \code{verbose=FALSE}, returns a list containing the following components:
#'   \describe{
#'     \item{\code{where}}{A description of how the marginal effects have been evaluated.}
#'     \item{\code{dEdq}}{The marginal effect. The partial derivative of the expected count with respect to a variable \code{q} in the \code{X} and or \code{Z} matrix, or the difference in expectation if \code{q} is binary.}
#'     \item{\code{se}}{The standard errors of the marginal effects evaluated numerically and using a Jacobian via the delta method.}
#'   }
#' 
#' @details 
#' The function computes marginal effects for zero-truncated Poisson or negative binomial regression models. 
#' It handles different model types; `oneinflmodel` for one-inflated models, and `truncmodel` for standard count models.
#' The marginal effects are evaluated at either all data points and averaged (`AE`, the default), at the sample means of the variables (`EM`), or at a custom case.
#' The marginal effects for dummy variables are actually the differences in expected outcomes for values of the dummy of 1 and 0. 
#' The marginal effects are displayed along with their statistical significance, evaluated based on the chosen `at` parameter.
#'
#' @examples
#' df <- data.frame(x = rnorm(100), z = rnorm(100), y = rpois(100, lambda = 5))
#' model <- oneinfl(y ~ x | z, df = df, dist = "Poisson")
#' margins(model, df, at = "AE") # Average Effect
#' margins(model, df, at = "EM", verbose=FALSE) # Effect at Means, suppress printing
#' margins(model, df, at = list(x = 1, z = 0)) # Custom case
#'
#' @seealso 
#' \code{\link{dEdq_nb}}, \code{\link{dEdq_nb_noinfl}}, \code{\link{dEdq_pois}}, 
#' \code{\link{dEdq_pois_noinfl}}, \code{\link{model.frame}}, 
#' \code{\link{model.matrix}}, \code{\link[stats]{numericDeriv}}
#' @export

margins <- function(model, df, at = "AE", verbose = TRUE) {
  
  # Function to create table
  create_table <- function(margins, se) {
    z_value <- margins / se
    p_value <- 2 * (1 - pnorm(abs(z_value)))
    
    tabl <- cbind(
      Marginal.effects = margins,
      Std.Error = se,
      z_value = z_value,
      p.value = p_value
    )
    
    significance <- sapply(p_value, get_significance)
    result_df <- data.frame(tabl, significance)
    colnames(result_df)[5] <- ""
    return(result_df)
  }
  
  # Function to determine significance symbols
  get_significance <- function(p_value) {
    if (p_value < 0.001) {
      return("***")
    } else if (p_value < 0.01) {
      return("**")
    } else if (p_value < 0.05) {
      return("*")
    } else if (p_value < 0.1) {
      return(".")
    } else {
      return("")
    }
  }
  
  q <- list()
  
  b <- model$beta
  if (inherits(model, "oneinflmodel")) { g <- model$gamma }
  if (model$dist == "negbin") { a <- model$alpha }
  
  names(b) <- substring(names(b), 2)
  if (inherits(model, "oneinflmodel")) { names(g) <- substring(names(g), 2) }
  
  formula <- model$formula
  cleandata <- makeXZy(formula, df)
  X <- cleandata$X
  if (inherits(model, "oneinflmodel")) { Z <- cleandata$Z }
  
  # Determine which variables are dummies
  is.dummy.X <- function(X) { length(unique(X)) == 2 }
  if (inherits(model, "oneinflmodel")) { 
    is.dummy.Z <- function(Z) { length(unique(Z)) == 2 }
  }
  
  dummies <- colnames(X)[apply(X, 2, is.dummy.X)]
  
  # Add dummies from Z if they aren't already in X
  if (inherits(model, "oneinflmodel")) {
    dummies <- c(
      colnames(X)[apply(X, 2, is.dummy.X)], 
      colnames(Z)[apply(Z, 2, is.dummy.Z)][
        !colnames(Z)[apply(Z, 2, is.dummy.Z)] %in% 
          colnames(X)[apply(X, 2, is.dummy.X)]
      ]
    )
  }
  
  if (is.list(at)) {
    q$where <- "Marginal effect evaluated at: "
    for (i in 1:(length(at) - 1)) {
      q$where <- paste(q$where, names(at[i]), " = ", unlist(at[i]), ", ", sep = "")
    }
    q$where <- paste(q$where, names(at[length(at)]), " = ", unlist(at[length(at)]), sep = "")
    if (!all(names(at) %in% names(df))) {
      stop("variable names in 'at' must match those in the data")
    }
    for (i in 1:length(at)) {
      if (names(at[i]) %in% colnames(X)) { X <- '[<-'(X, , names(at[i]), unlist(at[i])) }
      if (names(at[i]) %in% colnames(Z)) { Z <- '[<-'(Z, , names(at[i]), unlist(at[i])) }
    }
  } else if (at == "EM") {
    q$where <- "Marginal effect evaluated at the sample means of the data"
    X <- matrix((colSums(X) / nrow(X)), 1, ncol(X), dimnames = list(1, colnames(X)))
    Z <- matrix((colSums(Z) / nrow(Z)), 1, ncol(Z), dimnames = list(1, colnames(Z)))
  } else if (at == "AE") {
    q$where <- "Marginal effect averaged over all data points"
  } else {
    stop("'at' must be 'AE' (average effect), 'EM' (effect at means), or a list of representative cases in which to evaluate the effect")
  }
  
  if (model$dist == "Poisson") {
    if (inherits(model, "oneinflmodel")) {
      q$dEdq <- colMeans(dEdq_pois(b, g, X, Z, dummies, formula))
      J <- as.matrix(colMeans(attr(numericDeriv(quote(dEdq_pois(b, g, X, Z, dummies, formula)), c("b", "g")), "gradient")))
      q$se <- sqrt(diag(J %*% model$vc %*% t(J)))
    } else if (inherits(model, "truncmodel")) {
      q$dEdq <- colMeans(dEdq_pois_noinfl(b, X, dummies, formula))
      J <- as.matrix(colMeans(attr(numericDeriv(quote(dEdq_pois_noinfl(b, X, dummies, formula)), "b"), "gradient")))
      q$se <- sqrt(diag(J %*% model$vc %*% t(J)))
    }
  }
  
  if (model$dist == "negbin") {
    if (inherits(model, "oneinflmodel")) {
      q$dEdq <- colMeans(dEdq_nb(b, g, a, X, Z, dummies, formula))
      J <- as.matrix(colMeans(attr(numericDeriv(quote(dEdq_nb(b, g, a, X, Z, dummies, formula)), c("b", "g", "a")), "gradient")))
      q$se <- sqrt(diag(J %*% model$vc %*% t(J)))
    } else if (inherits(model, "truncmodel")) {
      q$dEdq <- colMeans(dEdq_nb_noinfl(b, a, X, dummies, formula))
      J <- as.matrix(colMeans(attr(numericDeriv(quote(dEdq_nb_noinfl(b, a, X, dummies, formula)), c("b", "a")), "gradient")))
      q$se <- sqrt(diag(J %*% model$vc %*% t(J)))
    }
  }
  
  if (inherits(model, "oneinflmodel")) {
    names(q$se) <- names(q$dEdq) 
    bnames <- names(b)[-1]
    gnames <- names(g)[-1]
    gnames <- gnames[!gnames %in% bnames]
    q$dEdq <- q$dEdq[c(bnames, gnames)]
    q$se <- q$se[c(bnames, gnames)]
  } else if (inherits(model, "truncmodel")) {
    names(q$se) <- names(q$dEdq) 
    bnames <- names(b)[-1]
    q$dEdq <- q$dEdq[c(bnames)]
    q$se <- q$se[c(bnames)]
  }
  
  # Creating table
  margins_table <- create_table(q$dEdq, q$se)
  
  # Conditional printing of results
  if (verbose) {
    cat("Call:\n")
    cat(paste("formula: ", deparse(model$formula), "\n"))
    cat(paste("distribution: ", model$dist, "\n"))
    
    cat("\nMarginal effects:\n")
    print(margins_table, digits = 4)
    
    cat(paste("\nSignif. codes:  0 ***' 0.001 **' 0.01 *' 0.05 .' 0.1  ' 1\n"))
  }
  
  if (verbose == FALSE) {return(q)} # Return a list instead of printing
}