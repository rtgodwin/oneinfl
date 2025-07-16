#' Predicted Expected Response for One-Inflated or Truncated Models
#'
#' Calculates the predicted expected response for a model fitted using 
#' \code{\link{oneinfl}} or \code{\link{truncreg}}.
#'
#' @param object An object of class `oneinflmodel`
#' @param ... Additional argument `df`, a data frame used to calculate the expected
#' value of the response variable.
#'
#' @return
#' A numeric vector of predicted expected responses for the observations in \code{df}.
#'
#' @details
#' This function computes the expected response based on the fitted model. The computation 
#' differs depending on the distribution. For \code{Poisson (OIPP)}, predicted values are
#' computed using \code{\link{E_pois}}. For \code{Negative Binomial (OIZTNB)}, predicted
#' values are computed using \code{\link{E_negbin}}.
#' 
#' @seealso
#' \code{\link{oneinfl}} for fitting one-inflated models.
#' \code{\link{E_pois}}, \code{\link{E_negbin}}, for the expected value calculations.
#'
#' @examples
#' # Example usage
#' df <- data.frame(x = rnorm(100), z = rnorm(100), y = rpois(100, lambda = 5))
#' model <- oneinfl(y ~ x | z, df = df, dist = "Poisson")
#' predict(model, df = df)
#'
#' @export

predict.oneinflmodel <- function(object, ...) {
  args <- list(...)
  df <- args$df

  b <- object$beta
  g <- object$gamma
  if (object$dist == "negbin") { a <- object$alpha }
  formula <- object$formula
  cleandata <- makeXZy(formula, df)
  X <- cleandata$X
  Z <- cleandata$Z
  
  if (object$dist == "negbin") {
    return(E_negbin(b, g, a, X, Z))
  }
  if (object$dist == "Poisson") {
    return(E_pois(b, g, X, Z))
  }
}