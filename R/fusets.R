#' @title Univariate time series fusion using regresssion based prediction
#' 
#' @description Univariate time series fusion using regresssion based prediction
#' 
#' @param xts original x time series; object of class "ts"
#' @param yts original y time series 
#' @param weight vector of weights (one for each correlation pair)
#' @param max_ewf maximum ewf to be considered)
#' @param steps steps at which the ewf is tested
#' @param plot optional plot the weight optimization plot. default=FALSE (no plot) 
#'
#' @author Johannes Reiche
#' 
#' @export 

fusets<- function(xts,yts,xi,yi,weight){
  
  y <- as.vector(na.omit(xi))
  x <- as.vector(na.omit(yi))
  #lm model
  fit=lm(y~poly(x,1),weight=weight)
  
  #prediction
  new_data <- na.omit(as.zoo(yts))
  
  pre <- predict(fit, data.frame(x=new_data))
  
  mts <- merge.zoo(xts, yts)
  mts[,1][!is.na(mts[,2])] <- as.zoo(pre)
  xfus <- mts[,1]

  return(xfus)
}