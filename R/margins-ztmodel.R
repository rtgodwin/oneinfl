margins.ztmodel <- function(model, data, at="AE") {
  q <- list()
  
  b <- model$beta
  if (model$dist == "negbin") {a <- model$alpha}
  
  names(b) <- substring(names(b), 2)
  
  formula <- model$formula
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  
  # determine which variables are dummies
  is.dummy <- function(X) {length(unique(X)) == 2}
  dummies <- colnames(data[-1])[apply(data[-1], 2, is.dummy)]
  
  if (is.list(at)) {
    q$where <- "Marginal effect evaluated at: "
    for(i in 1:(length(at) - 1)) {
      q$where <- paste(q$where, names(at[i]), " = ", unlist(at[i]), ", ", sep = "")
    }
    q$where <- paste(q$where, names(at[length(at)]), " = ", unlist(at[length(at)]), sep = "")
    if (!all(names(at) %in% names(data))) {stop("variable names in 'at' must match those in the data")}
    for(i in 1:length(at)) {
      if(names(at[i]) %in% colnames(X)) {X <- '[<-'(X, , names(at[i]), unlist(at[i]))}
    }
  } else if (at == "EM") {
    q$where <- "Marginal effect evaluated at the sample means of the data"
    X <- matrix((colSums(X) / nrow(X)), 1, ncol(X), dimnames = list(1, colnames(X)))
  } else if (at == "AE") {q$where <- "Marginal effect averaged over all data points"
  } else {
    stop("'at' must be 'AE' (average effect), 'EM' (effect at means), or a list of representative cases in which to evaluate the effect")
  }
  
  if (model$dist == "Poisson") {
    q$dEdq <- colMeans(dEdq_pois_noinfl(b, X, dummies))
    J <- as.matrix(colMeans(attr(numericDeriv(quote(dEdq_pois_noinfl(b, X, dummies)), "b"), "gradient")))
    q$se <- sqrt(diag(J %*% model$vc %*% t(J)))
  }
  
  if (model$dist == "negbin") {
    q$dEdq <- colMeans(dEdq_nb_noinfl(b, a, X, dummies))
    J <- as.matrix(colMeans(attr(numericDeriv(quote(dEdq_nb_noinfl(b, a, X, dummies)), c("b", "a")), "gradient")))
    q$se <- sqrt(diag(J %*% model$vc %*% t(J)))
  }
  q
}
