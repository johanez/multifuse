#' @title Univariate time series fusion, using regresssion based prediction
#' 
#' @description Univariate time series fusion, using regresssion based prediction. x time series (xts) gets fused with y time series (yts). The results is the original x time series with additional observations at the observation times of the y time series. A detailed description can be found in Reiche et al. 2015 (Section 3.1.3. "Univariate time series fusion")
#' 
#' @param xts original x time series
#' @param yts original y time series 
#' @param xi interpolated overlapping x time series. xi is calculated with \link{calcRegWeight}
#' @param yi interpolated overlapping y time series. yi is calcualted with \link{calcRegWeight} 
#' @param w vector of (optimized) multifuse regression weight (one for each correlation pair). w is calculated with \link{calcRegWeight}
#'
#' @return xfus fusted time series
#' 
#' @author Johannes Reiche
#' 
#' @examples
#' ## load example data
#' load("fiji.Rdata") 
#' xts <- ndvi_ts #Landsat NDVI example time series
#' yts <- hv_ts   #ALSO PALSAR HV example time series
#' 
#' ## calculate regression weight
#' wc <- calcRegWeight(xts,yts)
#' wc     #show time series data frame
#' wc$w   #show inverse normalised magnitude (multifuse regression weight)
#' 
#' ## optimize regression weight
#' opt <- optimizeRegWeight(wc$xi,wc$yi,wc$w,max_ewf=10,steps=0.1,order=1,plot=TRUE)
#  ewf <- opt$ewf.optimized   #optimized ewf
#' 
#' ## fuse time series, using optimized regression weight
#' xfus <- fusets(xts,yts,wc$xi,wc$yi,wc$weight)
#' plot2ts(xts,yts)  #plot original x and y time series
#' plot2ts(xfus,yts) #plot fused x time series with original y time series
#' 
#' 
#' @export 

fusets<- function(xts,yts,xi,yi,w){
  
  y <- as.vector(na.omit(xi))
  x <- as.vector(na.omit(yi))
  
  #lm model
  fit=lm(y~poly(x,1),weight=w)
  
  #prediction
  new_data <- na.omit(as.zoo(yts))
  pre <- predict(fit, data.frame(x=new_data))
  
  mts <- merge.zoo(xts, yts)
  mts[,1][!is.na(mts[,2])] <- as.zoo(pre)
  xfus <- mts[,1]

  return(xfus)
}