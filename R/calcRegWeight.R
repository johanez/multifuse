#' @title Calculate multifuse regression weight (for two time series)
#' 
#' @description Calculates the multifuse regression weight (inverse normalised magnitude) for two time series. A detailed description is given in Reiche et al. 2015 (Section 3.1.2 "Weighted time series correlation" and Section 3.1.2.1 "Deriviation of the time series regression weight").
#' @param x time series 1; object of class \link{ts}
#' @param y time series 2; object of class \link{ts}
#' 
#' @return A time series data frame with (1) overlapping x and y time series, (2) interpolated x and y time series (xi, yi), (3) x and y magnitude (xmag, ymag) and (4) multifuse regression weight (w). 
#' 
#' @author Johannes Reiche
#' 
#' @references Reiche, J., Verbesselt, J., Hoekman, D. H. & Herold, M. (2015): Fusing Landsat and SAR time series to detect deforestation in the tropics. Remote Sensing of Environment. 156, 276-293. DOI: 10.1016/j.rse.2014.10.001. \url{http://www.sciencedirect.com/science/article/pii/S0034425714003885} 
#' 
#' @examples
#' ## load example data
#' data(tsexample.rda)
#' xts <- ndvi_ts #Landsat NDVI example time series
#' yts <- hv_ts   #ALSO PALSAR HV example time series
#' 
#' ## calculate regression weight
#' wc <- calcRegWeight(xts,yts)
#' wc     #show time series data frame
#' wc$w   #show inverse normalised magnitude (multifuse regression weight)             
#' 
#' ## plot time series
#' plot2ts(xts,yts,lab_ts1="x",lab_ts2="y",title="Original x and y time series")
#' plot2ts(wc$x,wc$y,lab_ts1="x",lab_ts2="y",title="Overlapping x and y time series")   
#' plot2ts(wc$xi,wc$yi,lab_ts1="xi",lab_ts2="yi",title="Interpolated overlapping x and y time series") 
#' 
#' 
#' @import raster
#' @import tseries
#' @export 
#' 


calcRegWeight <-function(x,y) {
  
  calc.y_mag_for_x <- function(x,y){   
    #Calculating y magnitude for x time series observations (Reiche et al. 2015; Equation 3)
    
    tx <- time(x)[!is.na(x)]
    ty <- time(y)[!is.na(y)]
    xn <- na.remove(x)
    yn <- na.remove(y)
    
    #Y magnitudes for x time series observations
    Y_mag_xn <- tx[NA]
    #Y magnitudes for x time series observations (all steps)
    Y_mag_x <- x
    Y_mag_x[] <-NA
    
    for(i in 1:length(tx)){
      #Exception 1 (y-ts starts after xi, (xi-1:xi) this part is not overlapping)
      if(min(ty)>tx[i]) {
        #for x[i-1]:x[i] <- NA
        Y_mag_x[which(time(x) == tx[i])] <- NA
        #Exception 2 (y-ts ends before xi-1, (xi-1:xi) this part is not overlapping)
      } else if (max(ty)<tx[i]){
        Y_mag_x[which(time(x) == tx[i])] <- NA  
        #no exception
      } else{
        Y_mag_xn[i] <- yn[min(which(ty>tx[i]))] - yn[max(which(ty<tx[i]))]
        Y_mag_xn
        Y_mag_x[which(time(x) == tx[i])] <- Y_mag_xn[i]  
      }
    }
    return(Y_mag_x)
  }
  
  calc.w <- function(x,y,xmag,ymag){
    #Calculate regression weight (Reiche et al. 2015; Equation 4 - 8)
    
    calc.invertedmag <- function(w){
      #Calculate inverted normalised magnitude (Reiche et al. 2015; Equation 4 and 5)
      #Inverted magnitude (Reiche et al. 2015; Equation 4)
      w2 <- max(w[])-w[]+min(w[])
      #Normalised inverted magnitude (Reiche et al. 2015; Equation 5)
      w3 <- w2*(1/sum(w2))
      return(w3)
    }
    
    #get magnitudes for each observation
    xmag[!is.na(x)]<-NA
    ymag[!is.na(y)]<-NA
    #Calculate inverted normalised magnitude (Reiche et al. 2015; Equation 4 and 5)
    ixmag <- xmag
    iymag <- ymag
    ixmag[!is.na(ixmag)] <- calc.invertedmag(na.omit(xmag))
    iymag[!is.na(iymag)] <- calc.invertedmag(na.omit(ymag))
    #Normalise the inverted magnitude by number of x and y ts observations (Reiche et al. 2015; Equation 6 and 7)
    tl <- sum(length(na.omit(ixmag)),length(na.omit(iymag)))
    nxmag <- ixmag*length(na.omit(ixmag))/tl
    nymag <- iymag*length(na.omit(iymag))/tl
    #Combine time series (Reiche et al. 2015; Equation 8)
    w <- nxmag[]
    w[!is.na(nymag)]<-na.omit(nymag)
    w <- na.omit(w)    
    return(w)
  }
  
  #Define start and end of overlappping period  
  zx <- as.zoo(x)
  zy <- as.zoo(y)
  st <- max(start(na.omit(zx)),start(na.omit(zy)))
  en <- min(end(na.omit(zx)),end(na.omit(zy)))
  
  #Merge time series based on the longest time in order to extend the shorter time series to the length of the longer
  #This is ecessary because "window" does not exceed time series and this would result in two differently with varying length
  if(length(zx)>=length(zy)){
    zu <- merge(zx,zy)
    zx <- zu[,1]
    zy <- zu[,2]
  }else{
    zu <- merge(zy,zx)
    zx <- zu[,2]
    zy <- zu[,1]
  }
  
  #Subset ts
  zxw <- window(zx,start = st, end = en)
  zyw <- window(zy,start = st, end = en)
  
  #Interpolate ts and subset
  xi <- na.approx(x)
  yi <- na.approx(y)
  xiw <- window(xi,start = st, end = en)
  yiw <- window(yi,start = st, end = en)
  xiw[is.na(zxw)&is.na(zyw)] <- NA
  yiw[is.na(zxw)&is.na(zyw)] <- NA
  
  #Calculating y magnitude for x time series observations (Reiche et al. 2015; Equation 3)
  #Calculate global start and end 
  glob_st <- min(start(na.omit(zx)),start(na.omit(zy)))
  glob_en <- max(end(na.omit(zx)),end(na.omit(zy)))
  ye <- as.ts(window(y,start=glob_st,end=glob_en,extend=TRUE))
  #Calculate magnitude
  y_mag_for_x <- abs(calc.y_mag_for_x(x,ye))
  x_mag_for_y <- abs(calc.y_mag_for_x(ye,x))
  #Subset magnitude time series
  ymag <- window(y_mag_for_x,start = st, end = en)
  xmag <- window(x_mag_for_y,start = st, end = en)
  
  
  #Merge time series into a joint time series data frame
  dummy <- xiw
  zxw <- zxw[!is.na(dummy)] 
  zyw <- zyw[!is.na(dummy)] 
  xiw <- xiw[!is.na(dummy)] 
  yiw <- yiw[!is.na(dummy)] 
  xmag <- xmag[!is.na(dummy)] 
  ymag <- ymag[!is.na(dummy)] 
  
  wc <- merge(zxw,zyw,xiw,yiw,xmag,ymag,zxw[]<-NA)
  names(wc) <- c("x", "y", "xi", "yi","xmag", "ymag","w")  
  
  #Calculate regression weight (Reiche et al. 2015; Equation 4 - 8)
  wc$w <- calc.w(wc$x,wc$y,wc$xmag,wc$ymag)
  
  return(wc)
}