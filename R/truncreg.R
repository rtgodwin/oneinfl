#' Truncated Regression Model
#'
#' Fits a positive Poisson (PP) or zero-truncated negative binomial (ZTNB) regression model.
#'
#' @param formula A symbolic description of the model to be fitted.
#' @param df A data frame containing the variables in the model.
#' @param dist A character string specifying the distribution to use. Options are `"Poisson"` or `"negbin"`.
#' @param start Optional. A numeric vector of starting values for the optimization process. Defaults to `NULL`, in which case starting values are attempted to be chosen automatically.
#' @param method A character string specifying the optimization method to be passed to \code{\link[stats]{optim}}. Defaults to `"BFGS"`.
#' 
#' @return An object of class `"truncmodel"` containing the following components:
#'   \describe{
#'     \item{\code{beta}}{Estimated coefficients for the regression model.}
#'     \item{\code{alpha}}{Dispersion parameter (only for negative binomial distribution).}
#'     \item{\code{vc}}{Variance-covariance matrix of the estimated parameters.}
#'     \item{\code{logl}}{Log-likelihood of the fitted model.}
#'     \item{\code{dist}}{The distribution used for the model ("Poisson" or "negbin").}
#'     \item{\code{formula}}{The formula used for the model.}
#'   }
#' 
#' @details
#' This function fits a regression model for zero-truncated counts. Zero-truncated models are used when the count data does not include zeros, such as in cases where only positive counts are observed.
#' 
#' The function supports two distributions:
#' - `"Poisson"`: Zero-truncated Poisson regression.
#' - `"negbin"`: Zero-truncated negative binomial regression.
#' 
#' The function uses numerical optimization via \code{\link[stats]{optim}} to estimate the parameters.
#' 
#' @seealso
#' \code{\link{summary}} for summarizing the fitted model.
#' 
#' @examples
#' # Example usage
#' df <- data.frame(x = rnorm(100), y = rpois(100, lambda = 1) + 1)
#' model <- truncreg(y ~ x, df = df, dist = "Poisson")
#' summary(model)
#' 
#' @export

truncreg <- function(formula, df, dist = "negbin", start = NULL, method = "BFGS") {
  
  llpp <- function(param) {
    l <- as.vector(exp(X %*% param[1:kx]))
    if(max(y) > 170) {
      log.fac.y <- y * log(y) - y
      log.fac.y[y < 171] <- log(factorial(y)[y < 171])
    } else if (max(y) < 171) {log.fac.y <- log(factorial(y))}
    return(sum(y * log(l) - log(exp(l) - 1) - log.fac.y))
  }
  
  llztnb <- function(param) {
    li <- as.vector(exp(X %*% param[1:kx]))
    a  <- param[kx + 1]
    if(max(y) > 170) {
      log.fac.y <- y * log(y) - y
      log.fac.y[y < 171] <- log(factorial(y)[y < 171])
    } else if (max(y) < 171) {log.fac.y <- log(factorial(y))}
    ymax <- max(y)
    terms <- weights <- rep(0, ymax)
    for(ii in 1:ymax) {
      terms[ii] <- log(a + ii - 1)
    }
    weights[1] <- sum(y > 1)
    for(ii in 2:ymax) {
      weights[ii] <- sum(y > (ii - 1))
    }
    gterm <- sum(terms * weights)
    return(sum((y == 1) * log(a) - log.fac.y + a * log(a) + y * log(li) - (a + y) * log(a + li) - log(1 - (a / (a + li)) ^ a)) + gterm)
  }
  
  findstart <- function() {
    bs <- 2 / (kx * apply(X, 2, max))
    if(dist == "Poisson") {bs}
    else if (dist == "negbin") {c(bs, 0.5)}
  }
  
  z <- list()
  class(z) <- "truncmodel"
  z$formula <- formula
  
  cleandata <- makeXZy(formula, df)
  X <- cleandata$X
  y <- cleandata$y
  
  n <- length(y)
  kx <- NCOL(X)
  
  z$dist <- dist
  
  if(is.null(start)) {
    pstart <- findstart()
  } else {
    pstart <- start
  }
  
  if (dist == "Poisson") {
    fitp <- suppressWarnings(optim(fn=llpp, par=pstart, method=method, control=list(fnscale=-1, maxit=1000), hessian = T))
    if (fitp$convergence > 0)
      warning("optimization failed to converge")
    z$beta <- fitp$par[1:kx]
    z$vc <- -solve(as.matrix(fitp$hessian))
    colnames(z$vc) <- rownames(z$vc) <- paste("b", colnames(X), sep="")
    z$logl <- fitp$value
  } else if (dist == "negbin") {
    fitnb <- suppressWarnings(optim(fn = llztnb, par = pstart, method=method, control=list(fnscale=-1, maxit=5000), hessian = T))
    if (fitnb$convergence > 0)
      warning("optimization failed to converge")
    z$beta <- fitnb$par[1:kx]
    z$alpha <- as.numeric(fitnb$par[kx + 1])
    z$vc <- -solve(as.matrix(fitnb$hessian))
    colnames(z$vc) <- rownames(z$vc) <- c(paste("b", colnames(X), sep=""), "alpha")
    z$logl <- fitnb$value
  } else {stop("dist must be either Poisson or negbin")}
  names(z$beta) <- paste("b", colnames(X), sep="")
  z
}
