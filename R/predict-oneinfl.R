predict.oneinfl <- function(model, data, type = "response") {
  b <- model$beta
  if(class(model) == "oneinflmodel") {g <- model$gamma}
  if (model$dist == "negbin") {a <- model$alpha}
  formula <- model$formula
  cleandata <- makeXZy(formula, data)
  X <- cleandata$X
  if(class(model) == "oneinflmodel") {Z <- cleandata$Z}
  
  if(type == "response") {
    if(class(model) == "oneinflmodel") {
      if (model$dist == "negbin") {
        return(E_negbin(b, g, a, X, Z))
      }
      if (model$dist == "Poisson") {
        return(E_pois(b, g, X, Z))
      }
    }
    if(class(model) == "truncmodel") {
      if (model$dist == "negbin") {
        return(E_negbin_noinfl(b, a, X))
      }
      if (model$dist == "Poisson") {
        return(E_pois_noinfl(b, X))
      }
    }
  }
}