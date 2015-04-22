#' @title Optimize regression weight
#' 
#' @description Optimize multifuse regression weight. The optimized regression weight is determined by detecting the minimum p-value for ewf (exponential weight factor, weight^ewf) starting from ewf = 0. In case a local maximum p-value occurs, the ewf that results in the minimum p-value before the local maxima is selected. A detailed description can be found in Reiche et al. 2015 (Section 3.1.2.2 "Regression weight optimization")
#' 
#' @param xi interpolated x time series 
#' @param yi interpolated y time series 
#' @param w vector of multifuse regression weights (one for each correlation pair)
#' @param ewf_max maximum ewf to be considered. default=2.
#' @param ewf_steps steps at which the ewf is tested. The default=0.1.
#' @param plot optional plot the weight optimization plot. default=FALSE (no plot) 
#'
#' @return A list with (1) optimized ewf (ewf.optimized), (2) p-value found for optimized ewf (pvalue.optimized), and further paramter
#'
#' @author Johannes Reiche
#' 
#' @references Reiche, J., Verbesselt, J., Hoekman, D. H. & Herold, M. (2015): Fusing Landsat and SAR time series to detect deforestation in the tropics. Remote Sensing of Environment. 156, 276-293. DOI: 10.1016/j.rse.2014.10.001. \url{http://www.sciencedirect.com/science/article/pii/S0034425714003885} 
#' 
#' @examples
#' ## load example data
#' data(tsexample)
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
#' @import raster
#' @import tseries
#' @export 
#' 


optimizeRegWeight <- function(xi,yi,w,ewf_max=2,ewf_steps=0.1,order=1,plot=FALSE){
  
  vxi <- as.vector(xi)
  vyi <- as.vector(yi)
  
  #empty vectors for r-squared, p-value and ewf
  vrsquared <- rep(NA, 1.0)
  vpvalue <- rep(NA, 1.0)
  vewf <- rep(NA, 1.0)
  l <- 1
  
  #calculate r-squared, p-value for increasing exponent
  for(i in seq(from=0, to=ewf_max, by=ewf_steps)) {
    fit=lm(na.omit(vyi)~poly(na.omit(vxi),order),weight=w^i)
    vpvalue[l] <- anova(fit)$'Pr(>F)'[1]
    vrsquared[l] <- summary(fit)$r.squared
    vewf[l] <- i
    l <- l + 1
  }
  
  #bind vectors and derive extrema
  a <- cbind(vewf,vpvalue,vrsquared)
  pvalue.localmin <- a[which(diff(c(FALSE,diff(a[,2])>0,TRUE))>0)]
  pvalue.globalmin <- a[which(a[,2] == min(a[,2])),1]
  pvalue.local_peak <- a[which(diff(c(FALSE,diff(a[,2])>=0,FALSE))<0)]
  rsquared.globalmax <- a[which(a[,3] == max(a[,3])),1]
  
  #Determine optimized regression weight
  if(length(unique(vpvalue))==1) {
    ewf.optimized <- 0
  } else {
    if(length(pvalue.local_peak)==0){
      #in case no local pvalue peak is found -> optimized ewf that resulted in global minimum p-value
      ewf.optimized <- pvalue.globalmin
    } else {
      #in case local pvalue peak is found -> optimized ewf = ewf that resulted in minimum p-value before peak
      ewf.optimized <- pvalue.localmin[max(pvalue.localmin <= pvalue.local_peak[1])]
    }
  }
  
  #pvalue and rsquared at opitmized ewf
  pvalue.optimized <- formatC(vpvalue[which(vewf == ewf.optimized)],digits=10,format="f")
  rsquared.optimized <- formatC(vrsquared[which(vewf == ewf.optimized)],digits=10,format="f")
  
  #Create output list
  output <- list(
    ewf.optimized = ewf.optimized,
    pvalue.optimized = pvalue.optimized,  
    rsquared.optimized = rsquared.optimized,
    v.ewf = vewf,
    v.pvalue = vpvalue,
    v.rsquared = vrsquared,
    pvalue.localmin = pvalue.localmin,
    pvalue.globalmin = pvalue.globalmin,
    pvalue.local_peak = pvalue.local_peak,
    rsquared.globalmax = rsquared.globalmax
  )
  
  #Weight optimization plot
  if (plot==TRUE){
    p <- plot(vewf,vpvalue,axes=F, xlab=expression("exponential weight factor (ewf) of w"^"ewf"), ylab="",type="l",col="black",pch=".",cex=0.8)
    p <- p + points(vewf,vpvalue, pch = 19,cex=0.5)
    p <- axis(2, col="black",lwd=1)
    mtext(2,text="p-value",line=2)
    axis(1,vewf)
    par(new=T)
    p <- plot(vewf,vrsquared,axes=F, ylim=c(0,1),xlab="", ylab="",type="l",lty="dotted",col="blue",pch=".")
    p <- p + points(vewf,vrsquared, pch = 18,cex=0.8,col="blue")
    p <- axis(4, col="blue",lwd=1)
    mtext(4,text=expression('r'^2),line=2,col="blue",col.axis="blue")
    p <- axis(4, ylim=c(0,1),col="blue",col.axis="blue",lwd=1)
    abline(v = ewf.optimized,lty=3, col="red",lwd=2)
    title(paste("- Weight Optimization Plot -","\n","\n","optimized ewf=", formatC(ewf.optimized,digits=1,format="f")," (pvalue=",formatC(pvalue.optimized,digits=4)," rsquare: ", round(summary(fit)$r.squared,digits=4),")",sep=""),cex.main=1)
  }
  
  return(output)
}
