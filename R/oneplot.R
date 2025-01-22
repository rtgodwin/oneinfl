#' Plot Observed Data and Model Predictions
#'
#' Generates a bar plot of observed count data and overlays predicted values 
#' from one or more models fitted using \code{\link{oneinfl}} or \code{\link{truncreg}}.
#'
#' @param model1 The first fitted model object, either a one-inflated model (class `"oneinflmodel"`) 
#'   or a truncated model (class `"truncmodel"`).
#' @param model2 Optional. A second fitted model object, structured similarly to \code{model1}.
#' @param model3 Optional. A third fitted model object, structured similarly to \code{model1}.
#' @param model4 Optional. A fourth fitted model object, structured similarly to \code{model1}.
#' @param df A data frame containing the variables used in the models.
#' @param maxpred Optional. The maximum count value to include in the plot. Defaults to the maximum observed count.
#' @param ylimit Optional. The upper limit for the y-axis. Defaults to 1.1 times the highest observed frequency.
#' @param ccex Optional. A numeric value controlling the size of plot points and lines. Defaults to \code{1.5}.
#'
#' @return
#' A plot is generated but no values are returned.
#'
#' @details
#' This function visualizes observed count data as a bar plot and overlays predicted 
#' values from up to four models. The function automatically detects the type of model 
#' (Poisson or negative binomial; one-inflated or truncated) and adjusts the plot 
#' accordingly. Predictions are generated using the \code{\link{pred}} function.
#'
#' Model types are distinguished by different point and line styles:
#' \itemize{
#'   \item Poisson (PP): Dark magenta, triangle-down
#'   \item Zero-truncated negative binomial (ZTNB): Red, diamond
#'   \item One-inflated Poisson (OIPP): Green, triangle-up
#'   \item One-inflated zero-truncated negative binomial (OIZTNB): Blue, circle
#' }
#'
#' The legend in the top-right corner of the plot indicates the models displayed.
#'
#' @seealso
#' \code{\link{oneinfl}} for fitting one-inflated models.
#' \code{\link{truncreg}} for fitting truncated models.
#' \code{\link{pred}} for generating predictions used in the plot.
#'
#' @examples
#' # Example usage
#' df <- data.frame(x = rnorm(100), z = rnorm(100), y = rpois(100, lambda = 5) + 1)
#' model1 <- oneinfl(y ~ x | z, df = df, dist = "Poisson")
#' model2 <- truncreg(y ~ x, df = df, dist = "negbin")
#' oneplot(model1, model2, df = df, maxpred = 10)
#'
#' @export

oneplot <- function(model1, model2, model3, model4, df, maxpred, ylimit, ccex) {
  
  plotpp <- function(model, df, maxpred) {
    preds <- pred(model, df, maxpred)
    points(x = df.bar[,1], y = preds, pch = 25, col = "darkmagenta", cex = ccex)
    lines(x = df.bar[,1], y = preds, col = "darkmagenta", lwd = ccex)
    leg <<- c(leg, "PP")
    cols <<- c(cols, "darkmagenta")
    pchs <<- c(pchs, 25)
  }
  
  plotztnb <- function(model, df, maxpred) {
    preds <- pred(model, df, maxpred)
    points(x = df.bar[,1], y = preds, pch = 18, col = "red", cex = ccex)
    lines(x = df.bar[,1], y = preds, col = "red", lwd = ccex)
    leg <<- c(leg, "ZTNB")
    cols <<- c(cols, "red")
    pchs <<- c(pchs, 18)
  }
  
  plotoipp <- function(model, df, maxpred) {
    preds <- pred(model, df, maxpred)
    points(x = df.bar[,1], y = preds, pch = 17, col = "green", cex = ccex)
    lines(x = df.bar[,1], y = preds, col = "green", lwd = ccex)
    leg <<- c(leg, "OIPP")
    cols <<- c(cols, "green")
    pchs <<- c(pchs, 17)
  }
  
  plotoiztnb <- function(model, df, maxpred) {
    preds <- pred(model, df, maxpred)
    points(x = df.bar[,1], y = preds, pch = 16, col = "blue", cex = ccex)
    lines(x = df.bar[,1], y = preds, col = "blue", lwd = ccex)
    leg <<- c(leg, "OIZTNB")
    cols <<- c(cols, "blue")
    pchs <<- c(pchs, 16)
  }
  
  formula <- model1$formula
  cleandata <- makeXZy(formula, df)
  y <- cleandata$y
  
  if (missing(maxpred)) {
    maxpred <- max(y)
  }
  
  if (missing(ylimit)) {
    ylimit <- max(tabulate(y)) * 1.1
  }
  
  if (missing(ccex)) {
    ccex <- 1.5
  }
  
  df.bar <- barplot(tabulate(y)[1:maxpred], names = 1:maxpred, xlab = "count", ylab = "frequency", col = "gray", ylim = c(0, ylimit))
  leg <- "actual data"
  cols <- "gray"
  pchs <- 15
  
  if (inherits(model1, "truncmodel") & model1$dist == "Poisson") {
    plotpp(model1, df, maxpred)
  } else if (inherits(model1, "truncmodel") & model1$dist == "negbin") {
    plotztnb(model1, df, maxpred)
  } else if (inherits(model1, "oneinflmodel") & model1$dist == "Poisson") {
    plotoipp(model1, df, maxpred)
  } else if (inherits(model1, "oneinflmodel") & model1$dist == "negbin") {
    plotoiztnb(model1, df, maxpred)
  }
  
  if (!missing(model2)) {
    if (inherits(model2, "truncmodel") & model2$dist == "Poisson") {
      plotpp(model2, df, maxpred)
    } else if (inherits(model2, "truncmodel") & model2$dist == "negbin") {
      plotztnb(model2, df, maxpred)
    } else if (inherits(model2, "oneinflmodel") & model2$dist == "Poisson") {
      plotoipp(model2, df, maxpred)
    } else if (inherits(model2, "oneinflmodel") & model2$dist == "negbin") {
      plotoiztnb(model2, df, maxpred)
    }
  }
  
  if (!missing(model3)) {
    if (inherits(model3, "truncmodel") & model3$dist == "Poisson") {
      plotpp(model3, df, maxpred)
    } else if (inherits(model3, "truncmodel") & model3$dist == "negbin") {
      plotztnb(model3, df, maxpred)
    } else if (inherits(model3, "oneinflmodel") & model3$dist == "Poisson") {
      plotoipp(model3, df, maxpred)
    } else if (inherits(model3, "oneinflmodel") & model3$dist == "negbin") {
      plotoiztnb(model3, df, maxpred)
    }
  }
  
  if (!missing(model4)) {
    if (inherits(model4, "truncmodel") & model4$dist == "Poisson") {
      plotpp(model4, df, maxpred)
    } else if (inherits(model4, "truncmodel") & model4$dist == "negbin") {
      plotztnb(model4, df, maxpred)
    } else if (inherits(model4, "oneinflmodel") & model4$dist == "Poisson") {
      plotoipp(model4, df, maxpred)
    } else if (inherits(model4, "oneinflmodel") & model4$dist == "negbin") {
      plotoiztnb(model4, df, maxpred)
    }
  }
  
  legend("topright", legend = leg, col = cols, pch = pchs, cex = 1)
}