pred <- function(model, data, maxpred) {
  
  b <- model$beta
  
  formula <- model$formula
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  y <- cleandata$y
  
  if(missing(maxpred)) {
    maxpred = max(y)
  }
  
  l <- exp(X %*% b)
  
  if(class(model) == "oneinflmodel") {
    g <- model$gamma
    Z <- cleandata$Z
    t <- exp(-Z %*% g)
  }
  
  if (model$dist == "negbin") {
    a <- model$alpha
    th <- l / a
  }
  
  if(class(model) == "truncmodel" & model$dist == "Poisson") {
    #PP
    pred <- rep(0, maxpred)
    for(pp in 1:(maxpred)){
      pred[pp] <- sum((l ^ (pp)) / ((exp(l) - 1) * factorial(pp)))
    }
    pred
  } else if(class(model) == "truncmodel" & model$dist == "negbin") {
    #ZTNB
    pred <- rep(0, maxpred)
    for(pp in 1:(maxpred)){
      pred[pp] <- sum((gamma(a + pp) / gamma(a) / gamma(pp + 1)) * ((1 / (1 + th)) ^ a) * ((th / (1 + th)) ^ pp) * (1 / (1 - (1 + th) ^ (-a))))
    }
    pred
  } else if(class(model) == "oneinflmodel" & model$dist == "Poisson") {
    #OIPP
    pred <- rep(0, maxpred)
    w <- -l * (exp(l) - l - 1) ^ -1 + (1 + l * (exp(l) - l - 1) ^ -1) * (1 + t) ^ -1
    pred[1] <- sum(w + (1 - w) * l / (exp(l) - 1))
    for(pp in 2:(maxpred)){
      pred[pp] <- sum((1 - w) * (l ^ (pp)) / ((exp(l) - 1) * factorial(pp)))
    }
    pred
  } else if(class(model) == "oneinflmodel" & model$dist == "negbin") {
    #OIZTNB
    pred <- rep(0, maxpred)
    P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
    L <- -P1 / (1 - P1)
    w <- L + (1 - L) / (1 + t)
    pred[1] <- sum(w + (1 - w) * a * ((1 / (1 + th)) ^ a) * (th / (1 + th - (1 + th) ^ (1 - a))))
    for(pp in 2:(maxpred)){
      pred[pp] <- sum((1 - w) * (gamma(a + pp) / gamma(a) / gamma(pp + 1)) * ((1 / (1 + th)) ^ a) * ((th / (1 + th)) ^ pp) * (1 / (1 - (1 + th) ^ (-a))))
    }
    pred
  }
}
