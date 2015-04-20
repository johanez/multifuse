#' @title Calculate probability for deforestation
#' 
#' @description Calculate probability for deforestation
#' @description i refers to the i-th future signal that can take values between i = 0,...,n. 
#' @param bayts,minN=NULL,maxN=NULL,maxt=NULL,chi_ts1=0.5,chi_ts2=0.5,start=NULL,end=NULL)
#' @param bayts "bayts" time series data frame
#' @param minN Minimum number of future observations to confirm the change; Default=NULL -> no minN
#' @param maxN Maximum number of future observations to confirm the change; Default=NULL -> no maxN
#' @param maxt Maximum time (years) to confirm the change; Default=NULL -> no maxt
#' @param chi_ts1 Threshold of PChange at which the change is confirmed in case of the last observation is from for time series 2; Default=0.5
#' @param chi_ts2 Threshold of PChange at which the change is confirmed in case of the last observation is from for time series 2; Default=0.5
#' @param start Start date of monitoring period; Default=NULL -> start is start of bayts time series
#' @param end End date of monitoring period; Default=NULL -> end is end of bayts time series
#' @author Johannes Reiche
#' @return Updated "bayts" time series data frame with changes if detected
#' @export 
#' 

######## ts.correlation_A3 ##################################
#1.Interpolation between observations for X and Y ts
# 1a) linear interpolation (linear)
# 1b) use modelled seasonality ts to interpolation between points
#     this calles function ts.get_season_ts (therfore the model coeefficient and maxbreaks have to be defined))
# 2. itype: one-way:  Interpolation for all data point for which Xexisits
#    itype: tow-ways:  Interpolation for all data point for which X or Y exisits
# 3. Also returns with the correlation for each cell with maxdt
#Output for each cell:
#"x", "y" : original x and y ts values
#"xi", "yi": original x/y with corresponding interpolated x/y value
#"xmag", "ymag" : magnitude difference between two observations for each x and y 
#"magdiff" : manitude difference for yi/xi 
#"xtimediff"/"ytimediff": relocated x/y with partner with NN (time) difference
#"timediff" : time differnece
#"indextimediff"  
#Settings for Approach A,B,C
#Approach 1: maxdt = NULL -> xi,yi
#Appraoch 2: maxdt = maxdt, itype = "one-way" -> xtimediff,ytimediff
#Appraoch 3: maxdt = NULL -> xi,yi and timediff and magdiff for weights
#weights need to be inverted -> the larger the timediff the smaller the weight


calcRegWeight <-function(x,y) {
  
  get.y_mag_x <- function(x,y){
    #function that returns magnitude of y-ts for x-ts according to reich006
    #function expects 2 ts with same length 
    #time steps for obs in x and y ts
    tx <- time(x)[!is.na(x)]
    ty <- time(y)[!is.na(y)]
    #corresponding values for time steps of x and y ts
    xn <- na.remove(x)
    yn <- na.remove(y)
    #Y magnitudes for x ts steps
    Y_mag_xn <- tx[NA]
    #Y magnitudes for x ts (all steps)
    Y_mag_x <- x
    Y_mag_x[] <-NA
    #for 2nd obs in x ts until end of ts
    for(i in 1:length(tx)){
      #Exception 1 (y-ts starts after xi, (xi-1:xi) this part is not overlapping)
      if(min(ty)>tx[i]) {
        #for x[i-1]:x[i] <- NA
        Y_mag_x[which(time(x) == tx[i])] <- NA
        #Exception 2 (y-ts ends before xi-1, (xi-1:xi) this part is not overlapping)
      } else if (max(ty)<tx[i]){
        #for x[i-1]:x[i] <- NA
        Y_mag_x[which(time(x) == tx[i])] <- NA  
        #no exception
      } else{
        Y_mag_xn[i] <- yn[min(which(ty>tx[i]))] - yn[max(which(ty<tx[i]))]
        Y_mag_xn
        #for x[i-1]:x[i] <- Y_mag_xn[i] 
        Y_mag_x[which(time(x) == tx[i])] <- Y_mag_xn[i]  
      }
    }
    return(Y_mag_x)
  }
  
  get.inverseweight <- function(w){
    #invert 
    w2 <- max(w[])-w[]+min(w[])
    #normalise to sum of weights
    #w2*100/sum(w2)
    #normalise to 0 - 1
    w3 <- w2*(1/sum(w2))
    return(w3)
  }
  
  calc.w <- function(x,y,xmag,ymag){
        #get magnitudes
        xmag[!is.na(x)]<-NA
        ymag[!is.na(y)]<-NA
        #get inverese weight
        ixmag <- get.inverseweight(na.omit(xmag))
        iymag <- get.inverseweight(na.omit(ymag))
        #na.omit(ixmag) + na.omit(ixmag)
        xmag[!is.na(xmag)] <- ixmag
        ymag[!is.na(ymag)] <- iymag
        #normalise by number of  in x/y ts
        tl <- sum(length(na.omit(xmag)),length(na.omit(ymag)))
        wxmagn <- xmag*length(na.omit(xmag))/tl
        wymagn <- ymag*length(na.omit(ymag))/tl
        #combine
        wxymag <- wxmagn[]
        wxymag[!is.na(wymagn)]<-na.omit(wymagn)
        w <- na.omit(wxymag)    
    return(w)
  }
  
  #concert days to zoo
  zx <- as.zoo(x)
  zy <- as.zoo(y)
  
  #start and end of overlappping period
  st <- max(start(na.omit(zx)),start(na.omit(zy)))
  en <- min(end(na.omit(zx)),end(na.omit(zy)))

  #merge ts based on the longest ts to extend the shorter one to the length of the longer
  #necessary since the command "window" does not exceed ts and this would result in two differently long ts
  if(length(zx)>=length(zy)){
    zu <- merge(zx,zy)
    zx <- zu[,1]
    zy <- zu[,2]
  }else{
    zu <- merge(zy,zx)
    zx <- zu[,2]
    zy <- zu[,1]
  }
  
  #subset ts
  zxw <- window(zx,start = st, end = en)
  zyw <- window(zy,start = st, end = en)
  
  #interpolate ts
  xi <- na.approx(x)
  yi <- na.approx(y)

  #subset interpolated ts
  xiw <- window(xi,start = st, end = en)
  yiw <- window(yi,start = st, end = en)
  
  xiw[is.na(zxw)&is.na(zyw)] <- NA
  yiw[is.na(zxw)&is.na(zyw)] <- NA
  
  ###################################################################
  # calculate magnitude difference between two "real" ts observations
  # all ts points in between two "real" ts observations are assingned
  # with magnitude difference
  glob_st <- min(start(na.omit(zx)),start(na.omit(zy)))
  glob_en <- max(end(na.omit(zx)),end(na.omit(zy)))
  #ye <- as.ts(window(y,start=start(x),end=end(x),extend=TRUE))
  ye <- as.ts(window(y,start=glob_st,end=glob_en,extend=TRUE))
  
  #get magnitude ts
  y_mag_forx <- abs(get.y_mag_x(x,ye))
  x_mag_fory <- abs(get.y_mag_x(ye,x))
  
  #subset magnitude difference ts  
  ymag <- window(y_mag_forx,start = st, end = en)
  xmag <- window(x_mag_fory,start = st, end = en)
  
  
  #Merge TS
  dummy <- xiw
  zxw <- zxw[!is.na(dummy)] 
  zyw <- zyw[!is.na(dummy)] 
  xiw <- xiw[!is.na(dummy)] 
  yiw <- yiw[!is.na(dummy)] 
  xmag <- xmag[!is.na(dummy)] 
  ymag <- ymag[!is.na(dummy)] 
  #dummy[!is.na(wc$xi)&!is.na(wc$yi)] 
  
  
  wc <- merge(zxw,zyw,xiw,yiw,xmag,ymag,zxw[]<-NA)
  names(wc) <- c("x", "y", "xi", "yi","xmag", "ymag","w")  

  wc$w <- calc.w(wc$x,wc$y,wc$xmag,wc$ymag)


  return(wc)
}