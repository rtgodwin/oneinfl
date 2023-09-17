oneinfl <- function(formula, data, dist = "negbin", start = NULL, method = "BFGS") {
  
  lloipp <- function(param) {
    l <- as.vector(exp(X %*% param[1:kx]))
    t <- as.vector(Z %*% param[(kx + 1):(kx + kz)])
    L <- -l / (exp(l) - l - 1)
    w <- L + (1 - L) / (1 + exp(-t))
    return(sum(log(1 - w) + (y == 1) * log(w / (1 - w) + l / (exp(l) - 1)) + (1 - (y==1)) * (y * log(l) - log(exp(l) - 1) - log(factorial(y)))))
  }
  
  lloiztnb <- function(param) {
    li <- as.vector(exp(X %*% param[1:kx]))
    t <- as.vector(Z %*% param[(kx + 1):(kx + kz)])
    a  <- param[kx + kz + 1]
    th <- li / a
    P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
    L <- -P1 / (1 - P1)
    w <- L + (1 - L) / (1 + exp(-t))
    ymax <- max(y)
    terms <- weights <- rep(0, ymax)
    for(ii in 1:ymax) {
      terms[ii] <- log(a + ii - 1)
    }
    weights[1] <- sum(y > 1)
    for(ii in 2:ymax) {
      weights[ii] <- sum(y > (ii - 1))
    }
    gterm <- terms * weights
    return(sum(log(1 - w) + (y == 1) * (log(w / (1 - w) + a * ((a / (a + li)) ^ a) * (li / (a + li - a * (1 + (li / a)) ^ (1 - a))))) + (1 - (y == 1)) * (a * log(a) - log(factorial(y)) + y * log(li) - (a + y) * log(a + li) - log(1 - (a / (a + li)) ^ a))) + sum(gterm))
  }
  
  findstart <- function() {
    bs <- 2 / (kx * apply(X, 2, max))
    gs <- 0.5 / (kz * apply(Z, 2, max))
    if(dist == "Poisson") {c(bs, gs)}
    else if (dist == "negbin") {c(bs, gs, 0.5)}
  }
  
  z <- list()
  class(z) <- "oneinflmodel"
  z$formula <- formula
  
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  Z <- cleandata$Z
  y <- cleandata$y
  
  n <- length(y)
  kx <- NCOL(X)
  kz <- NCOL(Z)
  
  z$dist <- dist
  
  if(is.null(start)) {
    pstart <- findstart()
  } else {
    pstart <- start
  }
  
  if (dist == "Poisson") {
    fitp <- suppressWarnings(optim(fn=lloipp, par=pstart, method=method, control=list(fnscale=-1, maxit=1000), hessian = T))
    if (fitp$convergence > 0)
      warning("optimization failed to converge")
    z$beta <- fitp$par[1:kx]
    z$gamma <- fitp$par[(kx + 1):(kx + kz)]
    z$vc <- -solve(as.matrix(fitp$hessian))
    colnames(z$vc) <- rownames(z$vc) <- c(paste("b", colnames(X), sep=""), paste("g", colnames(Z), sep=""))
    z$logl <- fitp$value
  } else if (dist == "negbin") {
    fitnb <- suppressWarnings(optim(fn = lloiztnb, par = pstart, method=method, control=list(fnscale=-1, maxit=5000), hessian = T))
    if (fitnb$convergence > 0)
      warning("optimization failed to converge")
    z$beta <- fitnb$par[1:kx]
    np <- kx + kz + 1
    z$gamma <- fitnb$par[(kx + 1):(kx + kz)]
    z$alpha <- as.numeric(fitnb$par[np])
    z$vc <- -solve(as.matrix(fitnb$hessian))
    colnames(z$vc) <- rownames(z$vc) <- c(paste("b", colnames(X), sep=""), paste("g", colnames(Z), sep=""), "alpha")
    z$logl <- fitnb$value
  } else {stop("dist must be either Poisson or negbin")}
  names(z$beta) <- paste("b", colnames(X), sep="")
  names(z$gamma) <- paste("g", colnames(Z), sep="")
  z
}
