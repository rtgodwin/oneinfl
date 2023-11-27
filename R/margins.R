margins <- function(model, data, at="AE") {
  
  # Function to create table
  create_table <- function(margins, se) {
    z_value <- margins / se
    p_value <- 2 * (1 - pnorm(abs(z_value)))
    
    tabl <- cbind(
      Marginal.effects = margins,
      Std.Error = se,
      z_value = z_value,
      p.value = p_value
    )
    
    significance <- sapply(p_value, get_significance)
    result_df <- data.frame(tabl, significance)
    colnames(result_df)[5] <- ""
    return(result_df)
  }
  
  # Function to determine significance symbols
  get_significance <- function(p_value) {
    if (p_value < 0.001) {
      return("***")
    } else if (p_value < 0.01) {
      return("**")
    } else if (p_value < 0.05) {
      return("*")
    } else if (p_value < 0.1) {
      return(".")
    } else {
      return("")
    }
  }
  
  q <- list()
  
  b <- model$beta
  if(class(model) == "oneinflmodel") {g <- model$gamma}
  if (model$dist == "negbin") {a <- model$alpha}
  
  names(b) <- substring(names(b), 2)
  if(class(model) == "oneinflmodel") {names(g) <- substring(names(g), 2)}
  
  formula <- model$formula
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  if(class(model) == "oneinflmodel") {Z <- cleandata$Z}

  # determine which variables are dummies
  is.dummy.X <- function(X) {length(unique(X)) == 2}
  if(class(model) == "oneinflmodel") {is.dummy.Z <- function(Z) {length(unique(Z)) == 2}}
  dummies <- colnames(data[-1])[apply(data[-1], 2, is.dummy.X)]
  # add dummies from Z if they aren't already in X
  if(class(model) == "oneinflmodel") {
    dummies <- c(colnames(data[-1])[apply(data[-1], 2, is.dummy.X)], colnames(data[-1])[apply(data[-1], 2, is.dummy.Z)][!colnames(data[-1])[apply(data[-1], 2, is.dummy.Z)] %in% colnames(data[-1])[apply(data[-1], 2, is.dummy.X)]])
  }
  
  if (is.list(at)) {
    q$where <- "Marginal effect evaluated at: "
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
    q$where <- "Marginal effect evaluated at the sample means of the data"
    X <- matrix((colSums(X) / nrow(X)), 1, ncol(X), dimnames = list(1, colnames(X)))
    Z <- matrix((colSums(Z) / nrow(Z)), 1, ncol(Z), dimnames = list(1, colnames(Z)))
  } else if (at == "AE") {q$where <- "Marginal effect averaged over all data points"
  } else {
    stop("'at' must be 'AE' (average effect), 'EM' (effect at means), or a list of representative cases in which to evaluate the effect")
  }
  
  if (model$dist == "Poisson") {
    if(class(model) == "oneinflmodel") {
      q$dEdq <- colMeans(dEdq_pois(b, g, X, Z, dummies))
      J <- as.matrix(colMeans(attr(numericDeriv(quote(dEdq_pois(b, g, X, Z, dummies)), c("b", "g")), "gradient")))
      q$se <- sqrt(diag(J %*% model$vc %*% t(J)))
    } else if(class(model) == "truncmodel") {
      q$dEdq <- colMeans(dEdq_pois_noinfl(b, X, dummies))
      J <- as.matrix(colMeans(attr(numericDeriv(quote(dEdq_pois_noinfl(b, X, dummies)), "b"), "gradient")))
      q$se <- sqrt(diag(J %*% model$vc %*% t(J)))
    }
  }
  
  if (model$dist == "negbin") {
    if(class(model) == "oneinflmodel") {
      q$dEdq <- colMeans(dEdq_nb(b, g, a, X, Z, dummies))
      J <- as.matrix(colMeans(attr(numericDeriv(quote(dEdq_nb(b, g, a, X, Z, dummies)), c("b", "g", "a")), "gradient")))
      q$se <- sqrt(diag(J %*% model$vc %*% t(J)))
    } else if(class(model) == "truncmodel") {
      q$dEdq <- colMeans(dEdq_nb_noinfl(b, a, X, dummies))
      J <- as.matrix(colMeans(attr(numericDeriv(quote(dEdq_nb_noinfl(b, a, X, dummies)), c("b", "a")), "gradient")))
      q$se <- sqrt(diag(J %*% model$vc %*% t(J)))
    }
  }
  
  names(q$se) <- names(q$dEdq) 
  bnames <- names(b)[-1]
  q$dEdq <- q$dEdq[c(bnames)]
  q$se <- q$se[c(bnames)]
  if(class(model) == "oneinflmodel") {
    gnames <- names(g)[-1]
    gnames <- gnames[!gnames %in% bnames]
    q$dEdq <- q$dEdq[c(bnames, gnames)]
    q$se <- q$se[c(bnames, gnames)]
  }
  
  # Creating table
  margins_table <- create_table(q$dEdq, q$se)
  
  # Printing the results
  cat("Call:\n")
  cat(paste("formula: ", deparse(model$formula), "\n"))
  cat(paste("distribution: ", model$dist, "\n"))
  
  cat("\nMarginal effects:\n")
  print(margins_table, digits = 4)
  
  cat(paste("\nSignif. codes:  0 `***' 0.001 `**' 0.01 `*' 0.05 `.' 0.1 ` ' 1\n"))
}
