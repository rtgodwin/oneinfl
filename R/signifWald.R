signifWald <- function(model, varname) {
  q <- list()
  bname <- paste("b", varname, sep="")
  gname <- paste("g", varname, sep="")
  b <- model$beta[which(names(model$beta) %in% bname)]
  g <- model$gamma[which(names(model$gamma) %in% gname)]
  Vhat <- model$vc
  Vhat <- Vhat[which(rownames(Vhat) %in% c(bname, gname)), which(colnames(Vhat) %in% c(bname, gname))]
  q$W <- as.numeric(t(as.vector(c(b,g))) %*% solve(Vhat) %*% as.vector(c(b,g)))
  q$pval <- 1 - pchisq(q$W, 2)
  return(q)
}
