#' Predicted Expected Response for One-Inflated or Truncated Models
#'
#' Calculates the predicted expected response for a model fitted using 
#' \code{\link{oneinfl}} or \code{\link{truncreg}}.
#'
#' @param model A fitted model object of class \code{"oneinflmodel"} or \code{"truncmodel"}.
#' @param df A data frame containing the predictor variables used in the model.
#' @param type A character string specifying the type of prediction. Currently, only 
#'   \code{"response"} is supported, which calculates the expected value of the response variable.
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
#' predict(model, df = df, type = "response")
#'
#' @export

predict.oneinflmodel <- function(model, df, type = "response") {
  b <- model$beta
  g <- model$gamma
  if (model$dist == "negbin") { a <- model$alpha }
  formula <- model$formula
  cleandata <- makeXZy(formula, df)
  X <- cleandata$X
  Z <- cleandata$Z
  
  if (type == "response") {
      if (model$dist == "negbin") {
        return(E_negbin(b, g, a, X, Z))
      }
      if (model$dist == "Poisson") {
        return(E_pois(b, g, X, Z))
      }
  }
}