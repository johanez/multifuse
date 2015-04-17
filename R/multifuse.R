#' @title Optimzie regression weight
#' 
#' @description The optimized regression weight is determined by detecting the minimum p-value for ewf (exponential weight factor, weight^ewf) starting from ewf = 0. In cse a local maximum p-value occurs, the ewf that results in the minimum p-value before the local maxima is selected.
#' 
#' @param xi interpolated time series 
#' @param yi interpolated time series 
#' @param weight vector of weights (one for each correlation pair)
#' @param max_ewf maximum ewf to be considered)
#' @param steps steps at which the ewf is tested
#' @param plot optional plot the weight optimization plot. default=FALSE (no plot) 
#'
#' @author Johannes Reiche
#' 
#' @export 


multifuse <- function(xts,yts,ewf=1,optimize=FALSE,max_ewf=2,steps=1,order=1,alpha=1,plot=FALSE,wpoints=FALSE){
  #step 1: Weighted time series correlation
  #step 1a: calculate regression weight
  wc <- calcRegWeight(xts,yts)
  
  #step 1b: optimize regression weight (if choosen)
  if(optimize==TRUE){
    opt <- optimizeRegWeight(wc$xi,wc$yi,wc$w,max_ewf=max_ewf,steps=steps,order=order,plot=plot)
    ewf <- opt$ewf.optimized
  } 
  
  #plot correlation plot
  if(plot==TRUE){
    r <- wRsquared(wc$xi,wc$yi,weight=wc$w^ewf,order=order,plot=plot,wpoints=wpoints,xlab="xts",ylab="yts")
  } else {
    r <- wRsquared(wc$xi,wc$yi,weight=wc$w^ewf,order=order,plot=plot,wpoints=wpoints,xlab="xts",ylab="yts")   
  }
  
  if(alpha >= r$p.value){
  #step 2: univariate time series fusion, using regression based prediction
  xfus <- fusets(xts,yts,wc$xi,wc$yi,wc$w^ewf)
  } else {
    xfus <- NA
  }
  
  output <- list(
    xfus = xfus,
    ewf = ewf,
    p.value = r$p.value,
    r.squared = r$r.squared   
  )
  
  return(output)
}

