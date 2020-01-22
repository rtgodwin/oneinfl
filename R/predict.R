predplot <- function(model, data, maxpred = 6) {
  b <- model$beta
  g <- model$gamma
  if (model$dist == "negbin") {a <- model$alpha}

  formula <- model$formula
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  Z <- cleandata$Z
  y <- cleandata$y

  l <- exp(X %*% b)
  t <- exp(-Z %*% g)

  if (model$dist == "Poisson") {
    w <- -l * (exp(l) - l - 1) ^ -1 + (1 + l * (exp(l) - l - 1) ^ -1) * (1 + t) ^ -1
    predcontrol = predtreated = rep(0,maxpred)
    predcontrol[1] <- sum(l/(exp(l)-1))
    predtreated[1] <- sum(w + (1-w)*l/(exp(l)-1))
    for(pp in 2:(maxpred)){
      predcontrol[pp] <- sum((l^(pp))/((exp(l)-1)*factorial(pp)))
      predtreated[pp] <- sum((1-w)*(l^(pp))/((exp(l)-1)*factorial(pp)))
    }
  }

  if (model$dist == "negbin") {
    th <- l / a
    P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
    L <- -P1 / (1 - P1)
    w <- L + (1 - L) / (1 + t)
    predcontrol = predtreated = rep(0,maxpred)
    predcontrol[1] <- sum((gamma(a + 1) / gamma(a) / gamma(2)) * ((1 / (1 + th)) ^ a) * (th / (1 + th)) * (1 / (1 - (1 + th) ^ (-a))))
    predtreated[1] <- sum(w + (1 - w) * a * ((1 / (1 + th)) ^ a) * (th / (1 + th - (1 + th) ^ (1 - a))))
    for(pp in 2:(maxpred)){
      predcontrol[pp] <- sum((gamma(a + pp) / gamma(a) / gamma(pp + 1)) * ((1 / (1 + th)) ^ a) * ((th / (1 + th)) ^ pp) * (1 / (1 - (1 + th) ^ (-a))))
      predtreated[pp] <- sum((1 - w) * (gamma(a + pp) / gamma(a) / gamma(pp + 1)) * ((1 / (1 + th)) ^ a) * ((th / (1 + th)) ^ pp) * (1 / (1 - (1 + th) ^ (-a))))
    }
  }

  probs = t(matrix(c(predcontrol,predtreated),(maxpred),2))

  pdf("fee.pdf", width=4, height=3, pointsize = 12, bg = "white")
  par(mar=c(2,4,1,0.1),cex=.75,cex.axis=.75,cex.lab=0.75)
  df.bar <- barplot(probs, beside=TRUE,names=1:maxpred,xlab="count",ylab="frequency",col=c("gray","gray50"), ylim = c(0, tabulate(y)[1]*1.1))
  points(x = df.bar[2,], y = tabulate(y), pch=16, col="blue",cex=1)
  legend("topright", legend = c("actual data", "treatment distribution", "control (counterfactual) distribution"), col=c("blue", "gray50","gray"), pch = c(16,15,15), cex = 0.8)

  dev.off()
}
