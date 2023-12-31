oneLRT <- function(mod0, mod1) {
  k0 <- nrow(mod0$vc)
  k1 <- nrow(mod1$vc)
  if (k1 > k0) {
    LRTstat <- -2 * (mod0$logl - mod1$logl)
    pval <- 1 - pchisq(LRTstat, (k1 - k0))
  } else if (k0 > k1) {
    LRTstat <- -2 * (mod1$logl - mod0$logl)
    pval <- 1 - pchisq(LRTstat, (k0 - k1))
  } else {
    stop("models must differ in their number of parameters")
  }
  list(LRTstat = LRTstat, pval = pval)
}
