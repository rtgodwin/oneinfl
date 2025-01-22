#' Likelihood Ratio Test for Nested Models
#'
#' Performs a likelihood ratio test (LRT) to compare two nested models estimated by 
#' \code{\link{oneinfl}} or \code{\link{truncreg}}. It calculates the LRT statistic 
#' and its associated p-value, testing whether the more complex model provides 
#' a significantly better fit to the data than the simpler model.
#'
#' @param mod0 A model object (typically the simpler model) estimated using 
#'   \code{\link{oneinfl}} or \code{\link{truncreg}}.
#' @param mod1 A model object (typically the more complex model) estimated using 
#'   \code{\link{oneinfl}} or \code{\link{truncreg}}.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{LRTstat}}{The likelihood ratio test statistic.}
#'     \item{\code{pval}}{The p-value associated with the test statistic, based 
#'       on a chi-squared distribution.}
#'   }
#'
#' @details
#' The function extracts the log-likelihoods and number of parameters from the two 
#' models. It then calculates the LRT statistic:
#' \deqn{-2 \cdot (\ell_0 - \ell_1)}
#' where \eqn{\ell_0} and \eqn{\ell_1} are the log-likelihoods of the simpler and 
#' more complex models, respectively. The degrees of freedom for the test are 
#' equal to the difference in the number of parameters between the models.
#'
#' The likelihood ratio test is commonly used to test for:
#' - \emph{Overdispersion}: Comparing a Poisson model to a negative binomial model.
#' - \emph{One-inflation}: Comparing a one-inflated model to a non-one-inflated model.
#'
#' @seealso
#' \code{\link{oneinfl}} for fitting one-inflated models.
#' \code{\link{truncreg}} for fitting zero-truncated models.
#' \code{\link[stats]{pchisq}} for the chi-squared distribution.
#'
#' @examples
#' 
#' # Example: One-inflation test
#' df <- data.frame(y = rpois(100, lambda = 5), x = rnorm(100), z = rnorm(100))
#' OIZTNB <- oneinfl(y ~ x | z, df = df, dist = "negbin")
#' ZTNB <- truncreg(y ~ x, df = df, dist = "negbin")
#' oneLRT(OIZTNB, ZTNB)
#'
#' # Example: Overdispersion test
#' OIPP <- oneinfl(y ~ x | z, df = df, dist = "Poisson")
#' oneLRT(OIZTNB, OIPP)
#' 
#' @export

oneLRT <- function(mod0, mod1) {
  k0 <- nrow(mod0$vc)
  k1 <- nrow(mod1$vc)
  if (k1 > k0) {
    LRTstat <- -2 * (mod0$logl - mod1$logl)
    pval <- 1 - pchisq(LRTstat, (k1 - k0))
  } else if (k0 > k1) {
    LRTstat <- -2 * (mod1$logl - mod0$logl)
    pval <- 1 - pchisq(LRTstat, (k0 - k1))
  } else {
    stop("models must differ in their number of parameters")
  }
  list(LRTstat = LRTstat, pval = pval)
}
