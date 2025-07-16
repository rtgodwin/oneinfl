#' Plot Observed Data and Model Predictions
#'
#' Generates a bar plot of observed count data and overlays predicted values 
#' from one or more models fitted using \code{\link{oneinfl}} or \code{\link{truncreg}}.
#'
#' @param model1 The first fitted model object, either a one-inflated model (class `"oneinflmodel"`) 
#'   or a truncated model (class `"truncmodel"`).
#' @param model2 Optional. A second fitted model object, structured similarly to \code{model1}.
#' @param model3 Optional. A third fitted model object.
#' @param model4 Optional. A fourth fitted model object.
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
#' (Poisson or negative binomial; one-inflated or truncated) and adjusts the plot accordingly. 
#' Predictions are generated using the \code{\link{pred}} function.
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
#' model2 <- truncreg(y ~ x, df = df, dist = "Poisson")
#' oneplot(model1, model2, df = df, maxpred = 10)
#' 
#' @export

oneplot <- function(model1, model2 = NULL, model3 = NULL, model4 = NULL, df, maxpred = NULL, ylimit = NULL, ccex = 1.5) {
  plot_model <- function(model, df, maxpred, col, pch, label, leg, cols, pchs) {
    preds <- pred(model, df, maxpred)
    points(x = df.bar[,1], y = preds, pch = pch, col = col, cex = ccex)
    lines(x = df.bar[,1], y = preds, col = col, lwd = ccex)
    leg <- c(leg, label)
    cols <- c(cols, col)
    pchs <- c(pchs, pch)
    list(leg = leg, cols = cols, pchs = pchs)
  }
  
  formula <- model1$formula
  cleandata <- makeXZy(formula, df)
  y <- cleandata$y
  
  if (is.null(maxpred)) maxpred <- max(y)
  if (is.null(ylimit)) ylimit <- max(tabulate(y)) * 1.1
  
  df.bar <- barplot(tabulate(y)[1:maxpred], names = 1:maxpred, xlab = "count", ylab = "frequency", col = "gray", ylim = c(0, ylimit))
  
  leg <- "actual data"
  cols <- "gray"
  pchs <- 15
  
  models <- list(model1, model2, model3, model4)
  model_labels <- c("PP", "ZTNB", "OIPP", "OIZTNB")
  model_colors <- c("darkmagenta", "red", "green", "blue")
  model_pchs <- c(25, 18, 17, 16)
  
  for (model in models) {
    if (!is.null(model)) {
      if (inherits(model, "truncmodel") & model$dist == "Poisson") {
        res <- plot_model(model, df, maxpred, model_colors[1], model_pchs[1], model_labels[1], leg, cols, pchs)
      } else if (inherits(model, "truncmodel") & model$dist == "negbin") {
        res <- plot_model(model, df, maxpred, model_colors[2], model_pchs[2], model_labels[2], leg, cols, pchs)
      } else if (inherits(model, "oneinflmodel") & model$dist == "Poisson") {
        res <- plot_model(model, df, maxpred, model_colors[3], model_pchs[3], model_labels[3], leg, cols, pchs)
      } else if (inherits(model, "oneinflmodel") & model$dist == "negbin") {
        res <- plot_model(model, df, maxpred, model_colors[4], model_pchs[4], model_labels[4], leg, cols, pchs)
      }
      leg <- res$leg
      cols <- res$cols
      pchs <- res$pchs
    }
  }
  
  legend("topright", legend = leg, col = cols, pch = pchs, cex = 1)
}