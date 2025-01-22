#' Generate Predicted Frequencies for Models
#'
#' Computes the predicted frequencies for a specified model, given the observed data and 
#' maximum count value. This function supports models estimated using \code{\link{oneinfl}} 
#' or \code{\link{truncreg}}.
#'
#' @param model A fitted model object, either a one-inflated model (class `"oneinflmodel"`) 
#'   or a truncated model (class `"truncmodel"`).
#' @param df A data frame containing the variables used in the model.
#' @param maxpred Optional. The maximum count value for which predictions are generated. 
#'   Defaults to the maximum observed count in \code{df}.
#'
#' @return
#' A numeric vector of predicted frequencies for count values from 1 to \code{maxpred}.
#'
#' @details
#' The function computes predicted frequencies based on the type of model and its 
#' distribution:
#' \itemize{
#'   \item \strong{PP (Poisson, zero-truncated):} Computes predictions for a 
#'     zero-truncated Poisson model.
#'   \item \strong{ZTNB (Negative Binomial, zero-truncated):} Computes predictions 
#'     for a zero-truncated negative binomial model.
#'   \item \strong{OIPP (Poisson, one-inflated):} Computes predictions for a 
#'     one-inflated Poisson model.
#'   \item \strong{OIZTNB (Negative Binomial, one-inflated, zero-truncated):} Computes 
#'     predictions for a one-inflated, zero-truncated negative binomial model.
#' }
#'
#' The predictions are generated for count values from 1 to \code{maxpred}. For one-inflated 
#' models, the predictions account for the one-inflation probabilities.
#'
#' This is an internal function primarily used by \code{\link{oneplot}} for visualization purposes.
#'
#' @seealso
#' \code{\link{oneplot}} for visualizing observed and predicted frequencies.
#' \code{\link{oneinfl}} for fitting one-inflated models.
#' \code{\link{truncreg}} for fitting zero-truncated models.
#'
#' @keywords internal

pred <- function(model, df, maxpred) {
  
  b <- model$beta
  
  formula <- model$formula
  cleandata <- makeXZy(formula, df)
  X <- cleandata$X
  y <- cleandata$y
  
  if(missing(maxpred)) {
    maxpred = max(y)
  }
  
  l <- exp(X %*% b)
  
  if(inherits(model, "oneinflmodel")) {
    g <- model$gamma
    Z <- cleandata$Z
    t <- exp(-Z %*% g)
  }
  
  if (model$dist == "negbin") {
    a <- model$alpha
    th <- l / a
  }
  
  if(inherits(model, "truncmodel") & model$dist == "Poisson") {
    # PP
    pred <- rep(0, maxpred)
    for(pp in 1:(maxpred)){
      pred[pp] <- sum((l ^ (pp)) / ((exp(l) - 1) * factorial(pp)))
    }
    pred
  } else if(inherits(model, "truncmodel") & model$dist == "negbin") {
    # ZTNB
    pred <- rep(0, maxpred)
    for(pp in 1:(maxpred)){
      pred[pp] <- sum((gamma(a + pp) / gamma(a) / gamma(pp + 1)) * ((1 / (1 + th)) ^ a) * ((th / (1 + th)) ^ pp) * (1 / (1 - (1 + th) ^ (-a))))
    }
    pred
  } else if(inherits(model, "oneinflmodel") & model$dist == "Poisson") {
    # OIPP
    pred <- rep(0, maxpred)
    w <- -l * (exp(l) - l - 1) ^ -1 + (1 + l * (exp(l) - l - 1) ^ -1) * (1 + t) ^ -1
    pred[1] <- sum(w + (1 - w) * l / (exp(l) - 1))
    for(pp in 2:(maxpred)){
      pred[pp] <- sum((1 - w) * (l ^ (pp)) / ((exp(l) - 1) * factorial(pp)))
    }
    pred
  } else if(inherits(model, "oneinflmodel") & model$dist == "negbin") {
    # OIZTNB
    pred <- rep(0, maxpred)
    P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
    L <- -P1 / (1 - P1)
    w <- L + (1 - L) / (1 + t)
    pred[1] <- sum(w + (1 - w) * a * ((1 / (1 + th)) ^ a) * (th / (1 + th - (1 + th) ^ (1 - a))))
    for(pp in 2:(maxpred)){
      pred[pp] <- sum((1 - w) * (gamma(a + pp) / gamma(a) / gamma(pp + 1)) * ((1 / (1 + th)) ^ a) * ((th / (1 + th)) ^ pp) * (1 / (1 - (1 + th) ^ (-a))))
    }
    pred
  }
}