makeXZy <- function(formula, data) {
  if (formula[[3]] == 1) {
    X <- Z <- model.matrix(~1, data.frame(intercept = rep(1, nrow(data))))
    mf <- model.frame(formula = formula, data = data)
    y <- model.response(mf, "numeric")
    list(X = X, Z = Z, y = y)
  } else if (formula[[3]][[3]] == 1) {
    ff <- formula
    formula[[3]][1] <- call("+")
    ffc <- . ~ .
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    mf <- model.frame(formula = formula, data = data)
    mtX <- terms(ffc, data = data)
    X <- model.matrix(mtX, mf)
    Z <- model.matrix(~1, data.frame(intercept = rep(1, nrow(data))))
    y <- model.response(mf, "numeric")
    list(X = X, Z = Z, y = y)
  } else {
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
    mf <- model.frame(formula = formula, data = data)
    mtX <- terms(ffc, data = data)
    X <- model.matrix(mtX, mf)
    mtZ <- terms(ffz, data = data)
    mtZ <- terms(update(mtZ, ~.), data = data)
    Z <- model.matrix(mtZ, mf)
    y <- model.response(mf, "numeric")
    list(X = X, Z = Z, y = y)
  }
}
