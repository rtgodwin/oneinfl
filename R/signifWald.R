#' Wald Test for Significance of a Predictor Variable
#'
#' Performs a Wald test to evaluate the joint significance of a predictor variable in both 
#' the rate and one-inflation components of a model.
#'
#' @param model A fitted model object of class \code{"oneinflmodel"}.
#' @param varname A character string specifying the name of the predictor variable to test.
#'
#' @return
#' A list with the following components:
#' \describe{
#'   \item{\code{W}}{The Wald test statistic.}
#'   \item{\code{pval}}{The p-value associated with the test statistic, based 
#'     on a chi-squared distribution with 2 degrees of freedom.}
#' }
#'
#' @details
#' This function tests the null hypothesis that the coefficients for the specified predictor 
#' variable are jointly equal to zero in both the rate (\code{beta}) and one-inflation 
#' (\code{gamma}) components of the model. The test statistic is calculated as:
#' \deqn{W = \mathbf{c}^\top V^{-1} \mathbf{c}}
#' where \eqn{\mathbf{c}} is the vector of coefficients for the predictor in the rate and 
#' one-inflation components, and \eqn{V} is their variance-covariance matrix. The p-value is 
#' computed using a chi-squared distribution with 2 degrees of freedom.
#'
#' @seealso
#' \code{\link{oneinfl}} for fitting one-inflated models.
#' \code{\link{oneWald}} for a general Wald test of one-inflation parameters.
#'
#' @examples
#' # Example usage
#' set.seed(123)
#' df <- data.frame(x = rnorm(100), z = rnorm(100), y = rpois(100, lambda = 5))
#' model <- oneinfl(y ~ x | z, df = df, dist = "Poisson")
#' result <- signifWald(model, varname = "x")
#' print(result$W)    # Wald test statistic
#' print(result$pval) # p-value
#'
#' @export

signifWald <- function(model, varname) {
  q <- list()
  bname <- paste("b", varname, sep = "")
  gname <- paste("g", varname, sep = "")
  b <- model$beta[which(names(model$beta) %in% bname)]
  g <- model$gamma[which(names(model$gamma) %in% gname)]
  Vhat <- model$vc
  Vhat <- Vhat[which(rownames(Vhat) %in% c(bname, gname)), which(colnames(Vhat) %in% c(bname, gname))]
  q$W <- as.numeric(t(as.vector(c(b, g))) %*% solve(Vhat) %*% as.vector(c(b, g)))
  q$pval <- 1 - pchisq(q$W, 2)
  return(q)
}