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

optimizeRegWeight <- function(xi,yi,weight,max_ewf=2,steps=1,order=1,plot=FALSE){
  
  vxi <- as.vector(xi)
  vyi <- as.vector(yi)
  
  #empty vectors for R-square, p-value and ewf
  vrsquared <- rep(NA, 1.0)
  vpvalue <- rep(NA, 1.0)
  vewf <- rep(NA, 1.0)
  l <- 1
  
  #calculate R-square, p-value for increasing exponent
  for(i in seq(from=0, to=max_ewf, by=steps)) {
    fit=lm(na.omit(vyi)~poly(na.omit(vxi),order),weight=weight^i)
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

  #pvalue at opitmized ewf
  pvalue.optimized <- formatC(vpvalue[which(vewf == ewf.optimized)],digits=10,format="f")

  #create output list
  output <- list(
    ewf.optimized = ewf.optimized,
    pvalue.optimized = pvalue.optimized,    
    v.ewf = vewf,
    v.pvalue = vpvalue,
    v.rsquared = vrsquared,
    pvalue.localmin = pvalue.localmin,
    pvalue.globalmin = pvalue.globalmin,
    pvalue.local_peak = pvalue.local_peak,
    rsquared.globalmax = rsquared.globalmax
  )
  
  #weight optimization plot
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
  title(paste("- Weight Optimization Plot -","\n","\n","optimized ewf=", formatC(ewf.optimized,digits=1,format="f")," (p-value=",formatC(pvalue.optimized,digits=4),")",sep=""),cex.main=1)
  }
  
  return(output)
}
