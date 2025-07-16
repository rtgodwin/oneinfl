#' Wald Test for One-Inflation
#'
#' Performs a Wald test to evaluate the significance of the one-inflation parameters 
#' in a model estimated using \code{\link{oneinfl}}.
#'
#' @param model A model object of class \code{"oneinflmodel"} estimated using 
#'   \code{\link{oneinfl}}. The model must include one-inflation parameters (\code{gamma}) 
#'   and a variance-covariance matrix (\code{vc}).
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{W}}{The Wald test statistic.}
#'     \item{\code{pval}}{The p-value associated with the test statistic, based 
#'       on a chi-squared distribution.}
#'   }
#'
#' @details
#' The Wald test evaluates the null hypothesis that all one-inflation parameters 
#' (\code{gamma}) are equal to zero, indicating no one-inflation. The test statistic 
#' is calculated as:
#' \deqn{W = \gamma^\top V^{-1} \gamma}
#' where \eqn{\gamma} is the vector of one-inflation parameters and \eqn{V} is their 
#' variance-covariance matrix. The p-value is computed using a chi-squared distribution 
#' with degrees of freedom equal to the length of \eqn{\gamma}.
#'
#' This test is commonly used to determine whether a one-inflated model provides 
#' a significantly better fit than a non-one-inflated counterpart.
#'
#' @seealso
#' \code{\link{oneinfl}} for fitting one-inflated models.
#' \code{\link{oneLRT}} for a likelihood ratio test of nested models.
#' \code{\link[stats]{pchisq}} for the chi-squared distribution.
#'
#' @examples
#' # Example usage
#' df <- data.frame(y = rpois(100, lambda = 5), x = rnorm(100), z = rnorm(100))
#' OIZTNB <- oneinfl(y ~ x | z, df = df, dist = "negbin")
#' oneWald(OIZTNB)
#'
#' @export

oneWald <- function(model) {
  q <- list()
  g <- model$gamma
  Vhat <- model$vc
  Vhat <- Vhat[which(rownames(Vhat) %in% names(g)), which(colnames(Vhat) %in% names(g))]
  q$W <- as.numeric(t(as.vector(g)) %*% solve(Vhat) %*% as.vector(g))
  q$pval <- 1 - pchisq(q$W, length(g))
  return(q)
}
