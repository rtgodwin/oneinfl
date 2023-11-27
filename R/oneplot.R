oneplot <- function(model1, model2, model3, model4, data, maxpred, ylimit, ccex) {
  
  plotpp <- function(model, data, maxpred) {
    preds <- pred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch=25, col="darkmagenta",cex=ccex)
    lines(x = df.bar[,1], y = preds, col="darkmagenta",lwd=ccex)
    leg <<- c(leg, "PP")
    cols <<- c(cols, "darkmagenta")
    pchs <<- c(pchs, 25)
  }
  
  plotztnb <- function(model, data, maxpred) {
    preds <- pred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch=18, col="red",cex=ccex)
    lines(x = df.bar[,1], y = preds, col="red",lwd=ccex)
    leg <<- c(leg, "ZTNB")
    cols <<- c(cols, "red")
    pchs <<- c(pchs, 18)
  }
  
  plotoipp <- function(model, data, maxpred) {
    preds <- pred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch=17, col="green",cex=ccex)
    lines(x = df.bar[,1], y = preds, col="green",lwd=ccex)
    leg <<- c(leg, "OIPP")
    cols <<- c(cols, "green")
    pchs <<- c(pchs, 17)
  }
  
  plotoiztnb <- function(model, data, maxpred) {
    preds <- pred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch=16, col="blue",cex=ccex)
    lines(x = df.bar[,1], y = preds, col="blue",lwd=ccex)
    leg <<- c(leg, "OIZTNB")
    cols <<- c(cols, "blue")
    pchs <<- c(pchs, 16)
  }
  
  formula <- model1$formula
  cleandata <- makeXZy(formula, data)
  y <- cleandata$y
  
  if(missing(maxpred)) {
    maxpred = max(y)
  }
  
  if(missing(ylimit)) {
    ylimit = max(y) * 1.1
  }
  
  if(missing(ccex)) {
    ccex = 1.5
  }
  
  df.bar <- barplot(tabulate(y)[1:maxpred], names=1:maxpred, xlab="count", ylab="frequency", col="gray", ylim = c(0, ylimit))
  leg <- "actual data"
  cols <- "gray"
  pchs <- 15
  
  if(class(model1) == "truncmodel" & model1$dist == "Poisson") {
    plotpp(model1, data, maxpred)
  } else if(class(model1) == "truncmodel" & model1$dist == "negbin") {
    plotztnb(model1, data, maxpred)
  } else if(class(model1) == "oneinflmodel" & model1$dist == "Poisson") {
    plotoipp(model1, data, maxpred)
  } else if(class(model1) == "oneinflmodel" & model1$dist == "negbin") {
    plotoiztnb(model1, data, maxpred)
  }
  
  if(!missing(model2)) {
    if(class(model2) == "truncmodel" & model2$dist == "Poisson") {
      plotpp(model2, data, maxpred)
    } else if(class(model2) == "truncmodel" & model2$dist == "negbin") {
      plotztnb(model2, data, maxpred)
    } else if(class(model2) == "oneinflmodel" & model2$dist == "Poisson") {
      plotoipp(model2, data, maxpred)
    } else if(class(model2) == "oneinflmodel" & model2$dist == "negbin") {
      plotoiztnb(model2, data, maxpred)
    }
  }
  
  if(!missing(model3)) {
    if(class(model3) == "truncmodel" & model3$dist == "Poisson") {
      plotpp(model3, data, maxpred)
    } else if(class(model3) == "truncmodel" & model3$dist == "negbin") {
      plotztnb(model3, data, maxpred)
    } else if(class(model3) == "oneinflmodel" & model3$dist == "Poisson") {
      plotoipp(model3, data, maxpred)
    } else if(class(model3) == "oneinflmodel" & model3$dist == "negbin") {
      plotoiztnb(model3, data, maxpred)
    }
  }
  
  if(!missing(model4)) {
    if(class(model4) == "truncmodel" & model4$dist == "Poisson") {
      plotpp(model4, data, maxpred)
    } else if(class(model4) == "truncmodel" & model4$dist == "negbin") {
      plotztnb(model4, data, maxpred)
    } else if(class(model4) == "oneinflmodel" & model4$dist == "Poisson") {
      plotoipp(model4, data, maxpred)
    } else if(class(model4) == "oneinflmodel" & model4$dist == "negbin") {
      plotoiztnb(model4, data, maxpred)
    }
  }
  
  legend("topright", legend=leg, col=cols, pch=pchs, cex = 1)
}
