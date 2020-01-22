fee <- function(model, data, at="AE") {
  q <- list()

  b <- model$beta
  g <- model$gamma
  if (model$dist == "negbin") {a <- model$alpha}

  formula <- model$formula
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  Z <- cleandata$Z

  if (is.list(at)) {
    q$where <- "First exposure effect evaluated at: "
    for(i in 1:(length(at) - 1)) {
      q$where <- paste(q$where, names(at[i]), " = ", unlist(at[i]), ", ", sep = "")
    }
    q$where <- paste(q$where, names(at[length(at)]), " = ", unlist(at[length(at)]), sep = "")
    if (!all(names(at) %in% names(data))) {stop("variable names in 'at' must match those in the data")}
    for(i in 1:length(at)) {
      if(names(at[i]) %in% colnames(X)) {X <- '[<-'(X, , names(at[i]), unlist(at[i]))}
      if(names(at[i]) %in% colnames(Z)) {Z <- '[<-'(Z, , names(at[i]), unlist(at[i]))}
    }
  } else if (at == "EM") {
    q$where <- "First exposure effect evaluated at the sample means of the data"
    X <- colSums(X) / nrow(X)
    Z <- colSums(Z) / nrow(Z)
  } else if (at == "AE") {q$where <- "First exposure effect averaged over all data points"
  } else {
    stop("'at' must be 'AE' (average effect), 'EM' (effect at means), or a list of representative cases in which to evaluate the effect")
  }

  if (model$dist == "Poisson") {
    q$fee <- mean(fee_pois(b, g, X, Z))
    J <- as.matrix(colMeans(attr(numericDeriv(quote(fee_pois(b, g, X, Z)), c("b", "g")), "gradient")))
    q$sefee <- sqrt(t(J) %*% model$vc %*% J)
    q$treatment_visits <- sum(fee_pois(b, g, X, Z))

  } else if (model$dist == "negbin") {
    q$fee <- mean(fee_negbin(b, g, a, X, Z))
    J <- as.matrix(colMeans(attr(numericDeriv(quote(fee_negbin(b, g, a, X, Z)), c("b", "g", "a")), "gradient")))
    q$sefee <- sqrt(t(J) %*% model$vc %*% J)
  }
  q
}
