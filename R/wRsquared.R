#' @title Calculate weighted R-squared
#' 
#' @description Calculate weighted R-squared from two input time series, using weighted least squares. lm {stats} is used for fitting. 
#' 
#' @param x time series or vector (independent variable)
#' @param y time series or vector (dependent variable)
#' @param weigth an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. If non-NULL, weighted least squares is used with weights weights (that is, minimizing sum(w*e^2)); otherwise ordinary least squares is used. See also ‘Details’,
#' @param order degree of polynomial with which the relationship is modelled. default=1 (linear regression)
#' 
#' @param plot optional plot the correlation plot. default=FALSE (no plot)
#' @param wpoints size of the correlation points represent the weight. default=FALSE
#' @param xlab name for x (plot paramter)
#' @param ylab name for y (plot paramter)
#' @param xlim x limits (x1, x2). default=NULL (plot paramter)
#' @param ylab y limits (x1, x2). default=NULL (plot paramter)
#' 
#' @return rsquared and pvalue. See lm {stats} for description
#' 
#' @author Johannes Reiche
#' 
#' @export 

wRsquared <- function(x,y,weight=NULL,order=1,plot=FALSE,wpoints=FALSE,xlab="",ylab="",ylim=NULL,xlim=NULL){
  
  vx <- as.vector(x)
  vy <- as.vector(y)
  fit=lm(na.omit(vy)~poly(na.omit(vx),order),weight=weight)
  
  if(plot==TRUE){
    if(wpoints==TRUE){
      dat <- data.frame(X = x, Y = y, w = as.vector(weight*(2/max(weight))))
      plot(Y~X, data = dat, cex = w,pch=19,col = rgb(0,0,0, weight),ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab)
      points(Y~X, data = dat, cex = w,pch=1)
    } else {
      dat <- data.frame(X = x, Y = y, w = as.vector(weight))
      plot(Y~X, data = dat, ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab)  
    }

    r <- range(na.omit(vx))
    xx <- seq(r[1],r[2], length.out=250)
    predict(fit, data.frame(vx=xx))
    lines(xx, predict(fit, data.frame(vx=xx)), col='blue')
    title(paste("- Correlation Plot -","\n","\n"," r²: ", round(summary(fit)$r.squared,digits=4),";   p-value: ",formatC(anova(fit)$'Pr(>F)'[1],digits=10,format="f"),"   [n=",length(na.omit(x)),"]",sep=""),cex.main=1)    
  } 
  
  return(list(rsquared = summary(fit)$r.squared,pvalue = anova(fit)$'Pr(>F)'[1]))  
}