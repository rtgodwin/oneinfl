% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signifWald.R
\name{signifWald}
\alias{signifWald}
\title{Wald Test for Significance of a Predictor Variable}
\usage{
signifWald(model, varname)
}
\arguments{
\item{model}{A fitted model object of class \code{"oneinflmodel"}.}

\item{varname}{A character string specifying the name of the predictor variable to test.}
}
\value{
A list with the following components:
\describe{
  \item{\code{W}}{The Wald test statistic.}
  \item{\code{pval}}{The p-value associated with the test statistic, based 
    on a chi-squared distribution with 2 degrees of freedom.}
}
}
\description{
Performs a Wald test to evaluate the joint significance of a predictor variable in both 
the rate and one-inflation components of a model.
}
\details{
This function tests the null hypothesis that the coefficients for the specified predictor 
variable are jointly equal to zero in both the rate (\code{beta}) and one-inflation 
(\code{gamma}) components of the model. The test statistic is calculated as:
\deqn{W = \mathbf{c}^\top V^{-1} \mathbf{c}}
where \eqn{\mathbf{c}} is the vector of coefficients for the predictor in the rate and 
one-inflation components, and \eqn{V} is their variance-covariance matrix. The p-value is 
computed using a chi-squared distribution with 2 degrees of freedom.
}
\examples{
# Example usage
set.seed(123)
df <- data.frame(x = rnorm(100), z = rnorm(100), y = rpois(100, lambda = 5))
model <- oneinfl(y ~ x | z, df = df, dist = "Poisson")
result <- signifWald(model, varname = "x")
print(result$W)    # Wald test statistic
print(result$pval) # p-value

}
\seealso{
\code{\link{oneinfl}} for fitting one-inflated models.
\code{\link{oneWald}} for a general Wald test of one-inflation parameters.
}
