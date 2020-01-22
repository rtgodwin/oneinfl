oneinfl <- function(formula, data, dist = "negbin", start = NULL, method = "BFGS", starts = 1000, sdstarts = 3) {

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
    bs <- matrix(rnorm(kx * starts, 0, sdstarts), starts, kx)
    gs <- matrix(rnorm(kz * starts, 0, sdstarts), starts, kz)
    if (dist == "Poisson") {
      loglmat <- sapply(1:starts, function(q) {lloipp(c(bs[q, ], gs[q, ]))})
    } else if (dist == "negbin") {
      als <- runif(starts)
      loglmat <- sapply(1:starts, function(q) {lloiztnb(c(bs[q, ], gs[q, ], als[q]))})
    }
    bs <- bs[!is.infinite(loglmat), ]
    gs <- gs[!is.infinite(loglmat), ]
    if (dist == "negbin") {als <- als[!is.infinite(loglmat)]}
    loglmat <- loglmat[!is.infinite(loglmat)]
    if (formula[[3]] == 1) {
      bs <- bs[!is.nan(loglmat)]
      gs <- gs[!is.nan(loglmat)]
    } else {
      bs <- bs[!is.nan(loglmat), ]
      gs <- gs[!is.nan(loglmat), ]
    }
    if (dist == "negbin") {als <- als[!is.nan(loglmat)]}
    loglmat <- loglmat[!is.nan(loglmat)]
    coords <- which(loglmat == max(loglmat), arr.ind = T)
    if(length(loglmat) < 1) {stop("Failed to find starting values. Try increasing 'starts' from 1000")}
    if (dist == "Poisson") {
      if (formula[[3]] == 1) {c(bs[coords], gs[coords])
      } else {c(bs[coords, ], gs[coords, ])}
    }
    else if (dist == "negbin") {
      if (formula[[3]] == 1) {c(bs[coords], gs[coords], als[coords])
      } else {c(bs[coords, ], gs[coords, ], als[coords])}
    }
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

  if (dist == "Poisson") {
    fitp <- suppressWarnings(optim(fn=lloipp, par=findstart(), method=method, control=list(fnscale=-1, maxit=1000), hessian = T))
    if (fitp$convergence > 0)
      warning("optimization failed to converge")
    z$beta <- fitp$par[1:kx]
    z$gamma <- fitp$par[(kx + 1):(kx + kz)]
    z$vc <- -solve(as.matrix(fitp$hessian))
    colnames(z$vc) <- rownames(z$vc) <- c(paste("b", colnames(X), sep = "_"), paste("g", colnames(Z), sep = "_"))
    z$logl <- fitp$value
  } else if (dist == "negbin") {
    fitnb <- suppressWarnings(optim(fn = lloiztnb, par = findstart(), method=method, control=list(fnscale=-1, maxit=5000), hessian = T))
    if (fitnb$convergence > 0)
      warning("optimization failed to converge")
    z$beta <- fitnb$par[1:kx]
    np <- kx + kz + 1
    z$gamma <- fitnb$par[(kx + 1):(kx + kz)]
    z$alpha <- as.numeric(fitnb$par[np])
    z$vc <- -solve(as.matrix(fitnb$hessian))
    colnames(z$vc) <- rownames(z$vc) <- c(paste("beta", colnames(X), sep = "_"), paste("gamma", colnames(Z), sep = "_"), "alpha")
    z$logl <- fitnb$value
  } else {stop("dist must be either Poisson or negbin")}
  names(z$beta) <- paste("b", colnames(X), sep="")
  names(z$gamma) <- paste("g", colnames(Z), sep="")
  z
}
