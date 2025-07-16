#' Prepare Design Matrices and Response Vector
#'
#' Processes a model formula and a data frame to generate design matrices (`X` and `Z`) 
#' and a response vector (`y`) for regression models, including support for complex 
#' formulas with `|` operators.
#'
#' @param formula A symbolic description of the model, where the left-hand side specifies 
#'   the response variable and the right-hand side specifies predictors. 
#'   Formulas can include a `|` operator to separate predictors for different components of a model.
#' @param df A data frame containing the variables specified in the formula.
#' 
#' @return
#' A list containing the following components:
#' \describe{
#'   \item{\code{X}}{A design matrix for the main predictors.}
#'   \item{\code{Z}}{A design matrix for additional predictors (e.g., for a secondary process in 
#'     a two-component model).}
#'   \item{\code{y}}{The response vector extracted from the formula.}
#' }
#' 
#' @details
#' This function processes the formula to extract and construct:
#' - `X`: The main design matrix.
#' - `Z`: A secondary design matrix (if a `|` operator is used in the formula, separating components).
#' - `y`: The response variable.
#' 
#' It handles cases where the formula specifies:
#' - Only the main component (e.g., \code{y ~ x1 + x2}).
#' - A secondary component using the `|` operator (e.g., \code{y ~ x1 + x2 | z1 + z2}).
#' 
#' @seealso
#' \code{\link{model.matrix}}, \code{\link{model.frame}}, \code{\link{model.response}}
#' 
#' @export

makeXZy <- function(formula, df) {
  if (formula[[3]] == 1) {
    X <- Z <- model.matrix(~1, data.frame(intercept = rep(1, nrow(df))))
    mf <- model.frame(formula = formula, data = df)
    y <- model.response(mf, "numeric")
    list(X = X, Z = Z, y = y)
  }
  testerror <- try(formula[[3]][[3]] == 1, silent = TRUE)
  if (!inherits(testerror, "try-error")) {
    if (formula[[3]][[3]] == 1) {
      ff <- formula
      formula[[3]][1] <- call("+")
      ffc <- . ~ .
      ffc[[2]] <- ff[[2]]
      ffc[[3]] <- ff[[3]][[2]]
      mf <- model.frame(formula = formula, data = df)
      mtX <- terms(ffc, data = df)
      X <- model.matrix(mtX, mf)
      Z <- model.matrix(~1, data.frame(intercept = rep(1, nrow(df))))
      y <- model.response(mf, "numeric")
      list(X = X, Z = Z, y = y)
      }
  } 
  if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|"))) {
    ff <- formula
    formula[[3]][1] <- call("+")
    ffc <- . ~ .
    ffz <- ~.
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    ffz[[3]] <- ff[[3]][[3]]
    ffz[[2]] <- NULL
  } else {
      ffz <- ffc <- ff <- formula
      ffz[[2]] <- NULL
  }
    mf <- model.frame(formula = formula, data = df)
    mtX <- terms(ffc, data = df)
    X <- model.matrix(mtX, mf)
    mtZ <- terms(ffz, data = df)
    mtZ <- terms(update(mtZ, ~.), data = df)
    Z <- model.matrix(mtZ, mf)
    y <- model.response(mf, "numeric")
    list(X = X, Z = Z, y = y)
  
}