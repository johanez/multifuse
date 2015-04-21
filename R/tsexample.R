#' @title Landsat NDVI and ALSO PALSAR HV example time series 
#' 
#' @description Landsat NDVI (2005 - 2012, 11 observations) and ALOS PALSAR HV (2007 - 2010, 10 observations) time series for a pinus caribea planatation under very cloud conditions. Forest is harvested between July - September 2008 followed by replantation.  Two univariate time series objects of the class \link{ts}. 
#' 
#' @usage load("tsexample.rda")
#' 
#' @source Reiche, J., Verbesselt, J., Hoekman, D. H. & Herold, M. (2015): Fusing Landsat and SAR time series to detect deforestation in the tropics. Remote Sensing of Environment. 156, 276-293. DOI: 10.1016/j.rse.2014.10.001. \url{http://www.sciencedirect.com/science/article/pii/S0034425714003885} 
#' 
#' @author Johannes Reiche
#' 
#' @examples
#' ## load example data
#' load("tsexample.rda")
#' 
#' ## plot Landsat NDVI (ndvi) and ALSO PALSAR HV (hv) example time series
#' plot2ts(ndvi,hv,lab_ts1="Landsat NDVI",lab_ts2="ALOS PALSAR HV [dB]")
#' 
#' @export 

tsexample <- function(){
}