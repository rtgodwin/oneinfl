oneWald <- function(model) {
  q <- list()
  g <- model$gamma
  Vhat <- model$vc
  Vhat <- Vhat[which(rownames(Vhat) %in% names(g)), which(colnames(Vhat) %in% names(g))]
  q$W <- as.numeric(t(as.vector(g)) %*% solve(Vhat) %*% as.vector(g))
  q$pval <- 1 - pchisq(q$W, length(g))
  return(q)
}
