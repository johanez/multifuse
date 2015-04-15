#' @title Calculate weighted R-squared
#' 
#' @description Calculate weighted R-squared from two input time series,using weighted least squares
#' @description lm {stats} is used for fitting
#' @description Optional correlation plot
#' 
#' @param x time series or vector (independent variable)
#' @param y time series or vector (dependent variable)
#' @param weight an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. If non-NULL, weighted least squares is used with weights weights (that is, minimizing sum(w*e^2)); otherwise ordinary least squares is used. See also ‘Details’,
#' @param order degree of polynomial with which the relationship is modelled. default=1 (linear regression)
#' 
#' @param plot optional plot the correlation plot. default=FALSE (no plot)
#' @param xlab name for x (plot paramter)
#' @param ylab name for y (plot paramter)
#' @param xlim x limits (x1, x2). default=NULL (plot paramter)
#' @param ylab y limits (x1, x2). default=NULL (plot paramter)
#' @return r.squared and p.value. See lm {stats} for description (plot paramter)
#' 
#' @author Johannes Reiche
#' 
#' @export 

wRsquared <- function(x,y,weight=NULL,order=1,plot=FALSE,xlab="",ylab="",ylim=NULL,xlim=NULL){
  
  vx <- as.vector(x)
  vy <- as.vector(y)
  fit=lm(na.omit(vy)~poly(na.omit(vx),order),weight=weight)
  
  if(plot==TRUE){
    p <- plot(x,y,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
    r <- range(na.omit(vx))
    xx <- seq(r[1],r[2], length.out=250)
    predict(fit, data.frame(vx=xx))
    lines(xx, predict(fit, data.frame(vx=xx)), col='blue')
    title(paste("r²: ", round(summary(fit)$r.squared,digits=4),";   p-value: ",formatC(anova(fit)$'Pr(>F)'[1],digits=10,format="f"),"   [n=",length(na.omit(x)),"]",sep=""),cex.main=1)
  }
  
  return(list(r.squared = summary(fit)$r.squared,p.value = anova(fit)$'Pr(>F)'[1]))  
}


