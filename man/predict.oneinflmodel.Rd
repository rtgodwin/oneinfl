% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict-oneinfl.R
\name{predict.oneinflmodel}
\alias{predict.oneinflmodel}
\title{Predicted Expected Response for One-Inflated or Truncated Models}
\usage{
\method{predict}{oneinflmodel}(model, df, type = "response")
}
\arguments{
\item{model}{A fitted model object of class \code{"oneinflmodel"} or \code{"truncmodel"}.}

\item{df}{A data frame containing the predictor variables used in the model.}

\item{type}{A character string specifying the type of prediction. Currently, only
\code{"response"} is supported, which calculates the expected value of the response variable.}
}
\value{
A numeric vector of predicted expected responses for the observations in \code{df}.
}
\description{
Calculates the predicted expected response for a model fitted using
\code{\link{oneinfl}} or \code{\link{truncreg}}.
}
\details{
This function computes the expected response based on the fitted model. The computation
differs depending on the distribution. For \code{Poisson (OIPP)}, predicted values are
computed using \code{\link{E_pois}}. For \code{Negative Binomial (OIZTNB)}, predicted
values are computed using \code{\link{E_negbin}}.
}
\examples{
# Example usage
df <- data.frame(x = rnorm(100), z = rnorm(100), y = rpois(100, lambda = 5))
model <- oneinfl(y ~ x | z, df = df, dist = "Poisson")
predict(model, df = df, type = "response")

}
\seealso{
\code{\link{oneinfl}} for fitting one-inflated models.
\code{\link{E_pois}}, \code{\link{E_negbin}}, for the expected value calculations.
}
