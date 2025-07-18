% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roiztnb.r
\name{roiztnb}
\alias{roiztnb}
\title{Generate Random Counts from a One-Inflated Zero-Truncated Negative Binomial Process}
\usage{
roiztnb(b, g, alpha, X, Z)
}
\arguments{
\item{b}{A numeric vector of coefficients for the rate component.}

\item{g}{A numeric vector of coefficients for the one-inflation component.}

\item{alpha}{A numeric value representing the dispersion parameter for the negative binomial distribution.}

\item{X}{A matrix or data frame of predictor variables for the rate component.}

\item{Z}{A matrix or data frame of predictor variables for the one-inflation component.}
}
\value{
A numeric vector of simulated count data.
}
\description{
Simulates count data from a one-inflated, zero-truncated negative binomial (OIZTNB) process 
using specified parameters for the rate, one-inflation, and dispersion components.
}
\details{
This function generates count data from a one-inflated, zero-truncated negative binomial process. 
The process combines:
\itemize{
  \item A negative binomial distribution for counts greater than one.
  \item A one-inflation component that adjusts the probability of observing a count of one.
}

The algorithm:
\enumerate{
  \item Calculates the rate parameter (\eqn{\lambda}) as \eqn{\exp(X \cdot b)}.
  \item Computes the one-inflation probabilities (\eqn{\omega}) based on \eqn{Z \cdot g}.
  \item Computes the negative binomial dispersion parameter (\eqn{\theta = \lambda / \alpha}).
  \item Simulates counts for each observation:
    \itemize{
      \item Draws a random number to determine whether the count is one.
      \item Iteratively calculates probabilities for higher counts until the random number is matched.
    }
}

This function is useful for generating synthetic data for testing or simulation studies involving 
one-inflated, zero-truncated negative binomial models.
}
\examples{
# Example usage
set.seed(123)
X <- matrix(rnorm(100), ncol = 2)
Z <- matrix(rnorm(100), ncol = 2)
b <- c(0.5, -0.2)
g <- c(1.0, 0.3)
alpha <- 1.5
simulated_data <- roiztnb(b, g, alpha, X, Z)
print(simulated_data)

}
\seealso{
\code{\link{oneinfl}} for fitting one-inflated models.
}
