#' One-Inflated Regression Model
#'
#' Fits a one-inflated positive Poisson (OIPP) or one-inflated zero-truncated negative binomial (OIZTNB) regression model.
#'
#' @param formula A symbolic description of the model to be fitted. Variables before the pipe `|` link to the usual Poisson rate parameter, after the pipe link to the one-inflation parameter.
#' @param df A data frame containing the variables in the model.
#' @param dist A character string specifying the distribution to use. Options are `"Poisson"` or `"negbin"`.
#' @param start Optional. A numeric vector of starting values for the optimization process. Defaults to `NULL`, in which case starting values are attempted to be chosen automatically.
#' @param method A character string specifying the optimization method to be passed to \code{\link[stats]{optim}}. Defaults to `"BFGS"`.
#' 
#' @return An object of class `"oneinflmodel"` containing the following components:
#'   \describe{
#'     \item{\code{beta}}{Estimated coefficients for the rate component of the model.}
#'     \item{\code{gamma}}{Estimated coefficients for the one-inflation component of the model.}
#'     \item{\code{alpha}}{Dispersion parameter (only for negative binomial distribution).}
#'     \item{\code{vc}}{Variance-covariance matrix of the estimated parameters.}
#'     \item{\code{logl}}{Log-likelihood of the fitted model.}
#'     \item{\code{avgw}}{Average one-inflation probability.}
#'     \item{\code{absw}}{Mean absolute one-inflation probability.}
#'     \item{\code{dist}}{The distribution used for the model ("Poisson" or "negbin").}
#'     \item{\code{formula}}{The formula used for the model.}
#'   }
#' 
#' @details
#' This function fits a regression model for one-inflated counts. One-inflated models are used when there are an excess number of ones, relative to a Poisson or negative binomial process.
#' 
#' The function supports two distributions:
#' - `"Poisson"`: One-inflated Poisson regression.
#' - `"negbin"`: One-inflated negative binomial regression.
#' 
#' The function uses numerical optimization via \code{\link[stats]{optim}} to estimate the parameters.
#' 
#' @seealso
#' \code{\link{summary}} for summarizing the fitted model.
#' \code{\link{margins}} for calculating the marginal effects of regressors.
#' \code{\link{oneWald}} to test for no one-inflation.
#' \code{\link{signifWald}} for testing the joint significance of a single regressor that appears before and after the pipe `|`.
#' \code{\link{oneplot}} for plotting actual and predicted counts.
#' \code{\link{predict}} for expected response/dependent variable at each observation.
#' \code{\link{truncreg}} for fitting positive Poisson (PP) and zero-truncated negative binomial (ZTNB) models.
#' \code{\link{oneLRT}} to test for no one-inflation or no overdispersion using a nested PP, OIPP, or ZTNB model.
#' 
#' @examples
#' # Example usage
#' df <- data.frame(x = rnorm(100), z = rnorm(100), y = rpois(100, lambda = 1) + 1)
#' model <- oneinfl(y ~ x | z, df = df, dist = "Poisson")
#' summary(model)
#' margins(model, df)
#' oneWald(model)
#' predict(model, df=df)
#' 
#' @export

oneinfl <- function(formula, df, dist = "negbin", start = NULL, method = "BFGS") {
  
  lloipp <- function(param) {
    l <- as.vector(exp(X %*% param[1:kx]))
    t <- as.vector(Z %*% param[(kx + 1):(kx + kz)])
    L <- -l / (exp(l) - l - 1)
    w <- L + (1 - L) / (1 + exp(-t))
    if(max(y) > 170) {
      log.fac.y <- y * log(y) - y
      log.fac.y[y < 171] <- log(factorial(y)[y < 171])
    } else if (max(y) < 171) {log.fac.y <- log(factorial(y))}
    return(sum(log(1 - w) + (y == 1) * log(w / (1 - w) + l / (exp(l) - 1)) + (1 - (y==1)) * (y * log(l) - log(exp(l) - 1) - log.fac.y)))
  }
  
  lloiztnb <- function(param) {
    li <- as.vector(exp(X %*% param[1:kx]))
    t <- as.vector(Z %*% param[(kx + 1):(kx + kz)])
    a  <- param[kx + kz + 1]
    th <- li / a
    P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
    L <- -P1 / (1 - P1)
    w <- L + (1 - L) / (1 + exp(-t))
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
    gterm <- terms * weights
    return(sum(log(1 - w) + (y == 1) * (log(w / (1 - w) + a * ((a / (a + li)) ^ a) * (li / (a + li - a * (1 + (li / a)) ^ (1 - a))))) + (1 - (y == 1)) * (a * log(a) - log.fac.y + y * log(li) - (a + y) * log(a + li) - log(1 - (a / (a + li)) ^ a))) + sum(gterm))
  }
  
  findstart <- function() {
    bs <- 2 / (kx * apply(X, 2, max))
    gs <- 0.5 / (kz * apply(Z, 2, max))
    if(dist == "Poisson") {c(bs, gs)}
    else if (dist == "negbin") {c(bs, gs, 0.5)}
  }

  z <- list()
  class(z) <- "oneinflmodel"
  z$formula <- formula
  
  cleandata <- makeXZy(formula, df)
  X <- cleandata$X
  Z <- cleandata$Z
  y <- cleandata$y
  
  n <- length(y)
  kx <- NCOL(X)
  kz <- NCOL(Z)
  
  z$dist <- dist
  
  if(is.null(start)) {
    pstart <- findstart()
  } else {
    pstart <- start
  }
  
  if (dist == "Poisson") {
    fitp <- suppressWarnings(optim(fn=lloipp, par=pstart, method=method, control=list(fnscale=-1, maxit=1000), hessian = T))
    if (fitp$convergence > 0)
      warning("optimization failed to converge")
    z$beta <- fitp$par[1:kx]
    z$gamma <- fitp$par[(kx + 1):(kx + kz)]
    z$vc <- -solve(as.matrix(fitp$hessian))
    colnames(z$vc) <- rownames(z$vc) <- c(paste("b", colnames(X), sep=""), paste("g", colnames(Z), sep=""))
    z$logl <- fitp$value
    # Calculate the average one-inflation
    l <- exp(X %*% z$beta)
    t <- Z %*% z$gamma
    L <- -l / (exp(l) - l - 1)
    w <- L + (1 - L) / (1 + exp(-t))
    z$avgw <- mean(w)
    z$absw <- mean(abs(w))
  } else if (dist == "negbin") {
    fitnb <- suppressWarnings(optim(fn = lloiztnb, par = pstart, method=method, control=list(fnscale=-1, maxit=5000), hessian = T))
    if (fitnb$convergence > 0)
      warning("optimization failed to converge")
    z$beta <- fitnb$par[1:kx]
    np <- kx + kz + 1
    z$gamma <- fitnb$par[(kx + 1):(kx + kz)]
    z$alpha <- as.numeric(fitnb$par[np])
    z$vc <- -solve(as.matrix(fitnb$hessian))
    colnames(z$vc) <- rownames(z$vc) <- c(paste("b", colnames(X), sep=""), paste("g", colnames(Z), sep=""), "alpha")
    z$logl <- fitnb$value
    # Calculate the averagge one-inflation
    li <- exp(X %*% z$beta)
    t <- Z %*% z$gamma
    a  <- z$alpha
    th <- li / a
    P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
    L <- -P1 / (1 - P1)
    w <- L + (1 - L) / (1 + exp(-t))
    z$avgw <- mean(w)
    z$absw <- mean(abs(w))
  } else {stop("dist must be either Poisson or negbin")}
  names(z$beta) <- paste("b", colnames(X), sep="")
  names(z$gamma) <- paste("g", colnames(Z), sep="")
  z
}