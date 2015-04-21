#' @title Multi-sensor Time series Fusion (MulTiFuse) of two univariate time series
#' 
#' @description Pixel-based Multi-sensor Time series Fusion (MulTiFuse) of two univariat time series. "multifuse" wraps the core functions of the MulTiFuse approach. First, a weighted timeseries correlation is performed (\link{calcRegWeight}). To maximise the statistical significance of the correlation and to take exceptional cases into account, correlation weight optimization is done (\link{optimizeRegWeight}) before the relationship of the two time series is modelled through a weighted regression analysis \link{wRsquared}. The optimized regression model is utilized in a second step for regression-based prediction of time series observation to fuse the SAR and optical time series (\link{regfuse}). The output is a fused time series consisting of the original x time series observations and additional x observations at the observation times of the y time series. A detailed description of the MulTiFuse approach can be found in Reiche et. al. 2015 (see reference below).
#' 
#' @param xts Univariate time series 1; object of class \link{ts}
#' @param yts Univariate time series 2; object of class \link{ts}
#' @param ewf Fixed exponential weight factor (ewf) in case no regression weight optimization is done. Default value is ewf=1 (simple weight). ewf=0 (no weight) and ewf=2 (squared weight)
#' @param optimize Regression weight optimization to maximise the statistical significance of the correlation and to take exceptional cases into account
#' @param ewf_max Maximum ewf up to which the regression weight optimization is considered. Default value is 2.
#' @param ewf_steps ewf_steps steps at which the ewf is tested. The default value is 0.1
#' @param order Degree of polynomial with which the relationship is modelled. default=1 (linear regression)
#' @param plot Optional weight optimization plot and regression plot. 
#' @param alpha Optional p-value threshold up to which fusion is performed. In case pvalue is > alpha no fusion is performed. Default value is 1.
#'
#' @return (1) fused time series, (2) ewf, (3) pvalue and (4) rsquared
#' 
#' @author Johannes Reiche
#' 
#' @references Reiche, J., Verbesselt, J., Hoekman, D. H. & Herold, M. (2015): Fusing Landsat and SAR time series to detect deforestation in the tropics. Remote Sensing of Environment. 156, 276-293. DOI: 10.1016/j.rse.2014.10.001. \url{http://www.sciencedirect.com/science/article/pii/S0034425714003885} 
#' 
#' @examples
#' ## load example data
#' data(tsexample.rda)
#' xts <- ndvi #Landsat NDVI example time series
#' yts <- hv   #ALSO PALSAR HV example time series
#' plot2ts(xts,yts,lab_ts1="Landsat NDVI",lab_ts2="ALOS PALSAR HV [dB]")
#' 
#' ## apply multifuse
#' xfus <- multifuse(xts,yts,plot=TRUE)
#' 
#' xfus[[1]] #fused time series
#' xfus[[2]] #ewf
#' xfus[[3]] #pvalue
#' xfus[[4]] #rsquared
#' 
#' ##plot fused time series
#' plot2ts(xfus[[1]],yts,lab_ts1="Fused Landsat NDVI",lab_ts2="ALOS PALSAR HV [dB]")
#' 
#' @import raster
#' @import tseries
#' @export 

multifuse <- function(xts,yts,ewf=1,optimize=TRUE,ewf_max=2,ewf_steps=0.1,order=1,alpha=1,plot=FALSE,wpoints=FALSE){
  #step 1: Weighted time series correlation
  #step 1a: Calculate regression weight
  wc <- calcRegWeight(xts,yts)
  
  #step 1b: optimize regression weight (if choosen)
  if(optimize==TRUE){
    opt <- optimizeRegWeight(wc$xi,wc$yi,wc$w,ewf_max=ewf_max,ewf_steps=ewf_steps,order=order,plot=plot)
    ewf <- opt$ewf.optimized
  } 
  
  #plot correlation plot
  if(plot==TRUE){
    r <- wRsquared(wc$xi,wc$yi,weight=wc$w^ewf,order=order,plot=plot,wpoints=wpoints,xlab="xts",ylab="yts")
  } else {
    r <- wRsquared(wc$xi,wc$yi,weight=wc$w^ewf,order=order,plot=plot,wpoints=wpoints,xlab="xts",ylab="yts")   
  }
  
  if(alpha >= r$pvalue){
  #step 2: univariate time series fusion, using regression based prediction
  xfus <- regfuse(xts,yts,wc$xi,wc$yi,wc$w^ewf)
  } else {
    xfus <- NA
  }

  output <- list(
    xfus = xfus,
    ewf = ewf,
    pvalue = r$pvalue,
    rsquared = r$rsquared   
  )
  
  return(output)
}

