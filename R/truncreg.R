truncreg <- function(formula, data, dist = "negbin", start = NULL, method = "BFGS") {
  
  llpp <- function(param) {
    l <- as.vector(exp(X %*% param[1:kx]))
    return(sum(y * log(l) - log(exp(l) - 1) - log(factorial(y))))
  }
  
  llztnb <- function(param) {
    li <- as.vector(exp(X %*% param[1:kx]))
    a  <- param[kx + 1]
    ymax <- max(y)
    terms <- weights <- rep(0, ymax)
    for(ii in 1:ymax) {
      terms[ii] <- log(a + ii - 1)
    }
    weights[1] <- sum(y > 1)
    for(ii in 2:ymax) {
      weights[ii] <- sum(y > (ii - 1))
    }
    gterm <- sum(terms * weights)
    return(sum((y == 1) * log(a) - log(factorial(y)) + a * log(a) + y * log(li) - (a + y) * log(a + li) - log(1 - (a / (a + li)) ^ a)) + gterm)
  }
  
  findstart <- function() {
    bs <- 2 / (kx * apply(X, 2, max))
    if(dist == "Poisson") {bs}
    else if (dist == "negbin") {c(bs, 0.5)}
  }
  
  z <- list()
  class(z) <- "truncmodel"
  z$formula <- formula
  
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  y <- cleandata$y
  
  n <- length(y)
  kx <- NCOL(X)
  
  z$dist <- dist
  
  if(is.null(start)) {
    pstart <- findstart()
  } else {
    pstart <- start
  }
  
  if (dist == "Poisson") {
    fitp <- suppressWarnings(optim(fn=llpp, par=pstart, method=method, control=list(fnscale=-1, maxit=1000), hessian = T))
    if (fitp$convergence > 0)
      warning("optimization failed to converge")
    z$beta <- fitp$par[1:kx]
    z$vc <- -solve(as.matrix(fitp$hessian))
    colnames(z$vc) <- rownames(z$vc) <- paste("b", colnames(X), sep="")
    z$logl <- fitp$value
  } else if (dist == "negbin") {
    fitnb <- suppressWarnings(optim(fn = llztnb, par = pstart, method=method, control=list(fnscale=-1, maxit=5000), hessian = T))
    if (fitnb$convergence > 0)
      warning("optimization failed to converge")
    z$beta <- fitnb$par[1:kx]
    z$alpha <- as.numeric(fitnb$par[kx + 1])
    z$vc <- -solve(as.matrix(fitnb$hessian))
    colnames(z$vc) <- rownames(z$vc) <- c(paste("b", colnames(X), sep=""), "alpha")
    z$logl <- fitnb$value
  } else {stop("dist must be either Poisson or negbin")}
  names(z$beta) <- paste("b", colnames(X), sep="")
  z
}
