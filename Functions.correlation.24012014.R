
### Tools
####### rmse ##########
rmse <- function(obs, pred) sqrt(mean((obs-pred)^2))
###### inverse weights #################
#inverse weights and normalise to 0 - 1
get.inverseweight <- function(w){
  #for inverte weights
  #invert 
  w2 <- max(w[])-w[]+min(w[])
  #normalise to sum of weights
  #w2*100/sum(w2)
  #normalise to 0 - 1
  w3 <- w2*(1/sum(w2))
  return(w3)
}
###### get magweight ##############
#to be documented!
#provides magnitude weight for all pixel pairs, the higher the magnitude the lower the weight
get.magweight <- function(x,y,xmag,ymag,itype=c("one-way","two-ways")){
  if (itype=="one-way"){
    ymag[!is.na(y)]<-NA
    #get inveresd weight and normalised weight (0-1)
    #the larger the magnitude the lower the weight
    iymag <- get.inverseweight(na.omit(ymag))
    w <- iymag
  }else 
    if (itype=="two-ways"){
      #get magnitudes
      xmag[!is.na(x)]<-NA
      ymag[!is.na(y)]<-NA
      #na.omit(xmag)
      #na.omit(ymag)
      #get inverese weight
      ixmag <- get.inverseweight(na.omit(xmag))
      iymag <- get.inverseweight(na.omit(ymag))
      #na.omit(ixmag) + na.omit(ixmag)
      xmag[!is.na(xmag)] <- ixmag
      ymag[!is.na(ymag)] <- iymag
      #na.omit(xmag)
      #na.omit(ymag)
      #normalise by number of  in x/y ts
      tl <- sum(length(na.omit(xmag)),length(na.omit(ymag)))
      wxmagn <- xmag*length(na.omit(xmag))/tl
      wymagn <- ymag*length(na.omit(ymag))/tl
      #na.omit(wxmagn)
      #na.omit(wymagn)
      #combine
      wxymag <- wxmagn[]
      #na.omit(wxymag)
      wxymag[!is.na(wymagn)]<-na.omit(wymagn)
      w <- na.omit(wxymag)
    } else stop("Not a correct interpolation type is selected ('one-way' or 'two-ways') ")
  
  return(w)
}

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
#function that finds global maxima and minima
which.peaks <- function(x,partial=TRUE,localMin=FALSE){
  if (localMin){
    if (partial){
      which(diff(c(FALSE,diff(x)>0,TRUE))>0)
    }else {
      which(diff(diff(x)>0)>0)+1
    }
  }else {
    if (partial){
      which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
    }else {
      which(diff(diff(x)>=0)<0)+1
    }
  }
}
#Optimization of wt

opt_wt.plot <- function(xi,yi,wmag,max_wt_exp=10,steps=1,order=1){
  vxi <- as.vector(xi)
  vyi <- as.vector(yi)
  #empty vectors for R-square, p-value and exponent
  vR <- rep(NA, 1.0)
  vpVal <- rep(NA, 1.0)
  vwt_exp <- rep(NA, 1.0)
  l <- 1
  #calculate R-square, p-value for increasing exponent
  for(i in seq(from=0, to=max_wt_exp, by=steps)) {
    fit=lm(na.omit(vyi)~poly(na.omit(vxi),order),weight=wmag^i)
    vpVal[l] <- anova(fit)$'Pr(>F)'[1]
    vR[l] <- summary(fit)$r.squared
    vwt_exp[l] <- i
    l <- l + 1
  }
  
  #bind vectors and derive local minima & maxima
  a <- cbind(vwt_exp,vpVal,vR)
  #r-square local minima
  localmin_pVal <- a[which.peaks(a[,2],localMin=TRUE),1]
  globalmin_pVal <- a[which(a[,2] == min(a[,2])),1]
  localmin_rsquare <- a[which.peaks(a[,3],localMin=TRUE),1]
  nonpartial_localmin_rsquare <- a[which.peaks(a[,3],localMin=TRUE,partial=FALSE),1]
  nonpartial_localmax_pVal <- a[which.peaks(a[,2],localMin=FALSE,partial=FALSE),1]
  globalmax_rsquare <- a[which(a[,3] == max(a[,3])),1]
  
  if(length(unique(vpVal))==1) {
    wt_exp <- 0
  } else {
    #derive best fit (lowest p-value) under the side condition of a meaningful correlation
    #meaninful correlation is always before local minima of r-square
    if(length(nonpartial_localmax_pVal)==0){
      #if no non-partial local minima of r-squared is found,
      #the global minima of pVal
      wt_exp <- globalmin_pVal
    } else {
      #if non-partial local mimima found for r-squared,
      #then: largest local minima of pVal that is <= smalles non partial local minima
      wt_exp <- localmin_pVal[max(localmin_pVal <= nonpartial_localmax_pVal[1])]
    }
  }
  #corresponding p-values to exponent that was found to be the best
  pVal <- formatC(vpVal[which(vwt_exp == wt_exp)],digits=10,format="f")
  #create output list
  rval <- list(
    vwt_exp = vwt_exp,
    vpVal = vpVal,
    vrsquare = vR,
    localmin_pVal = localmin_pVal,
    globalmin_pVal = globalmin_pVal,
    localmin_rsquare = localmin_rsquare,
    nonpartial_localmin_rsquare = nonpartial_localmin_rsquare,
    globalmax_rsquare = globalmax_rsquare,
    nonpartial_localmax_pVal = nonpartial_localmax_pVal,
    wt_exp = wt_exp,
    pVal = pVal    
  )
  
  #plots
  p <- plot(vwt_exp,vpVal,axes=F, xlab=expression("Exponential weight factor (ewf) of wt"["2d"]^"   ewf"), ylab="",type="l",col="black",pch=".",cex=0.8)
  p <- p + points(vwt_exp,vpVal, pch = 19,cex=0.5)
  p <- axis(2, col="black",lwd=1)
  mtext(2,text="p-value",line=2)
  axis(1,vwt_exp)
  #plot y
  #cat ("Press [enter] to plot 2nd time-series")
  #line <- readline()
  par(new=T)
  p <- plot(vwt_exp,vR,axes=F, ylim=c(0,1),xlab="", ylab="",type="l",col="blue",pch=".")
  p <- p + points(vwt_exp,vR, pch = 19,cex=0.5,col="blue")
  p <- axis(4, col="blue",lwd=1)
  mtext(4,text=expression('r'^2),line=2,col="blue")
  p <- p + axis(4, ylim=c(0,1),col="blue",lwd=1)
  abline(v = nonpartial_localmin_rsquare,lty=3, col="blue")
  abline(v = localmin_pVal,lty=3, col="black",lwd=2)
  title(paste("Minimum p-value at ewf: ", formatC(wt_exp,digits=1,format="f")," [p-value: ",formatC(pVal,digits=10,format="f"),"]",sep=""),cex.main=1)
  return(rval)
}


opt_wt <- function(xi,yi,wmag,max_wt_exp=10,steps=1,order=1){
  vxi <- as.vector(xi)
  vyi <- as.vector(yi)
  #empty vectors for R-square, p-value and exponent
  vR <- rep(NA, 1.0)
  vpVal <- rep(NA, 1.0)
  vwt_exp <- rep(NA, 1.0)
  l <- 1
  #calculate R-square, p-value for increasing exponent
  for(i in seq(from=0, to=max_wt_exp, by=steps)) {
    fit=lm(na.omit(vyi)~poly(na.omit(vxi),order),weight=wmag^i)
    vpVal[l] <- anova(fit)$'Pr(>F)'[1]
    vR[l] <- summary(fit)$r.squared
    vwt_exp[l] <- i
    l <- l + 1
  }
  
  #bind vectors and derive local minima & maxima
  a <- cbind(vwt_exp,vpVal,vR)
  #r-square local minima
  localmin_pVal <- a[which.peaks(a[,2],localMin=TRUE),1]
  globalmin_pVal <- a[which(a[,2] == min(a[,2])),1]
  localmin_rsquare <- a[which.peaks(a[,3],localMin=TRUE),1]
  nonpartial_localmin_rsquare <- a[which.peaks(a[,3],localMin=TRUE,partial=FALSE),1]
  nonpartial_localmax_pVal <- a[which.peaks(a[,2],localMin=FALSE,partial=FALSE),1]
  globalmax_rsquare <- a[which(a[,3] == max(a[,3])),1]
  
  #in case all pairs have the same weight pVal will be equal for all
  #wt_exp, then all are local min and maxima and an error occurs
  #set wt_exp to 0
  #
  if(length(unique(vpVal))==1) {
    wt_exp <- 0
  } else {
    if(length(nonpartial_localmax_pVal)==0){
      #if no non-partial local minima of r-squared is found,
      #the global minima of pVal
      wt_exp <- globalmin_pVal
    } else {
      #if non-partial local mimima found for r-squared,
      #then: largest local minima of pVal that is <= smalles non partial local minima
      wt_exp <- localmin_pVal[max(localmin_pVal <= nonpartial_localmax_pVal[1])]
    }
  }
  #corresponding p-values to exponent that was found to be the best
  pVal <- formatC(vpVal[which(vwt_exp == wt_exp)],digits=30,format="f")
  
  #derive best fit (lowest p-value) under the side condition of a meaningful correlation
  #meaninful correlation is always before local minima of r-square
  
  #create output list
  rval <- list(
    vwt_exp = vwt_exp,
    vpVal = vpVal,
    vrsquare = vR,
    localmin_pVal = localmin_pVal,
    globalmin_pVal = globalmin_pVal,
    localmin_rsquare = localmin_rsquare,
    nonpartial_localmin_rsquare = nonpartial_localmin_rsquare,
    globalmax_rsquare = globalmax_rsquare,
    nonpartial_localmax_pVal = nonpartial_localmax_pVal,
    wt_exp = wt_exp,
    pVal = pVal    
  )
  return(rval)
}

######## ts.get.season_ts ##############
#(i)  Model seasonality for entire time-series including breaks (separate the ts in segments) and 
#(ii) Return ts of seasonal model
# adapted from JV script
ts.get_season_ts <- function(vi,order,maxbreaks){
  vipp <- bfastpp(vi, order = order) # you can adjust the order
  bp.vi <- breakpoints(response ~ 1, data = vipp, breaks = maxbreaks)
  
  #breakfactor - Factor Coding of Segmentations [strucchange]
  #Generates a factor encoding the segmentation given by a set of breakpoints.
  #vipp$seg <- breakfactor(bp.vi, breaks = nrofbreaks)
  vipp$seg <- breakfactor(bp.vi)
  
  #lm - Fitting Linear Models [stats]
  #watch the lm formula
  if (length(which(vipp$seg=="segment2"))==0){
    fm1 <- lm(response ~ (trend + harmon), 
              data = vipp) 
  } else {
    #lm - Fitting Linear Models [stats]
    fm1 <- lm(response ~ seg/(trend + harmon)-1, 
              data = vipp) 
  }
  
  zvi <- na.omit(as.zoo(vi))
  zvifit <- zoo(predict(fm1, vipp),time(zvi))
  
  ## Time of breakpoint
  bps <- breakpoints(bp.vi) ## two breakpoints
  abline(v = time(zvi)[bps$breakpoints], col = "red", lty = 2)
  
  ## now we can use the model for "interpolation"
  ## time steps at which we do the interpolation
  tts <- time(vi)
  newdat <- bfastpp(tts, order = order)
  
  ## set the segments manually
  ## labels and breaks has to have the same number otherwise it does not work
  ## number of breaks is set by the user before
  # wonder if this is dynaimc for different change mechanisms
  newdat$seg <- cut(newdat$time,
                    breaks = c(newdat$time[1], 
                               vipp$time[bp.vi$breakpoints], 
                               newdat$time[length(newdat$time)]),
                    labels = levels(vipp$seg),
                    include.lowest = TRUE)
  
  ## predict using the newdat object with predefined segments
  predict  <- predict(fm1, newdat)
  znew <- zoo(predict, time(tts))
  
  #return ts of seasonal model
  return(znew)
}

ts.get_season_ts.plot <- function(vi,order,maxbreaks){
  vipp <- bfastpp(vi, order = order) # you can adjust the order
  bp.vi <- breakpoints(response ~ 1, data = vipp, breaks = maxbreaks)
  
  #breakfactor - Factor Coding of Segmentations [strucchange]
  #Generates a factor encoding the segmentation given by a set of breakpoints.
  #vipp$seg <- breakfactor(bp.vi, breaks = nrofbreaks)
  vipp$seg <- breakfactor(bp.vi)
  #vipp$seg
  
  #lm - Fitting Linear Models [stats]
  #watch the lm formula
  if (length(which(vipp$seg=="segment2"))==0){
    fm1 <- lm(response ~ (trend + harmon), 
              data = vipp) 
  } else {
    #lm - Fitting Linear Models [stats]
    fm1 <- lm(response ~ seg/(trend + harmon)-1, 
              data = vipp) 
  }
  
  ## plotting
  zvi <- na.omit(as.zoo(vi))
  plot(zvi, col = "lightgrey", lwd = 3, ylab = "NDVI")
  zvifit <- zoo(predict(fm1, vipp),time(zvi))
  lines(zvifit, col = "lightblue", lwd = 2)
  
  ## Time of breakpoint
  bps <- breakpoints(bp.vi) ## two breakpoints
  abline(v = time(zvi)[bps$breakpoints], col = "red", lty = 2)
  
  ## now we can use the model for "interpolation"
  ## time steps at which we do the interpolation
  tts <- time(vi)
  newdat <- bfastpp(tts, order = order)
  
  ## set the segments manually
  ## labels and breaks has to have the same number otherwise it does not work
  ## number of breaks is set by the user before
  # wonder if this is dynaimc for different change mechanisms
  newdat$seg <- cut(newdat$time,
                    breaks = c(newdat$time[1], 
                               vipp$time[bp.vi$breakpoints], 
                               newdat$time[length(newdat$time)]),
                    labels = levels(vipp$seg),
                    include.lowest = TRUE)
  newdat$time[1]
  vipp$time[bp.vi$breakpoints]
  newdat$time[length(newdat$time)]
  levels(vipp$seg)
  
  ## predict using the newdat object with predefined segments
  predict  <- predict(fm1, newdat)
  znew <- zoo(predict, time(tts))
  #plot(znew, type = "p", cex = 0.1, col = "lightgrey", pch = 19)
  #points(vi, pch = 19, col = "red", cex = 0.3)
  
  #plot ts seasonal model and linear model fit 
  #(linear interpolation of actual datapoints for seasonal model)
  plot(znew,lwd = 3)
  points(zvi,pch = 19, col = "black", cex = 0.3)
  lines(zvi, col = "red", lwd = 2)
  abline(v = time(zvi)[bps$breakpoints], col = "red", lty = 2)
  #return ts of seasonal model
  return(znew)
}

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

ts.correlation_A3 <-function(x,y,maxdt=NULL,itype=c("one-way","two-ways"),interpolation_x=c("linear","season"),interpolation_y=c("linear","season"),order_x=2,order_y=2,maxbreaks_x=2,maxbreaks_y=2) {
  
  #Function that returns interpolated TS
  get.ts.interpolated <- function(x,interpolation=c("linear","season"),order=2,maxbreaks=2) {
    if(interpolation=="linear"){
      xi <- na.approx(x)
    } else 
      if (interpolation_x=="season"){
        #get season ts  
        xi <- ts.get_season_ts(x,order=order_x,maxbreaks=maxbreaks_x)
        #substitute points in season ts with real values
        xi[!is.na(as.zoo(x))] <- x[!is.na(as.zoo(x))]
      } else stop("Not a correct interpolation model is selected ('linear' or 'season') ")
    return(xi)
  }
  
  #concert days to zoo
  zx <- as.zoo(x)
  zy <- as.zoo(y)
  #start and end of overlappping period
  st <- max(start(na.omit(zx)),start(na.omit(zy)))
  en <- min(end(na.omit(zx)),end(na.omit(zy)))
  #include one valid observation exceeding the overlapping ts by maxdt days (max and min) 
  if(!is.null(maxdt)){
    zmaxdt=maxdt/365
    st <- st - zmaxdt
    en <- en + zmaxdt
  }
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
  xi <- get.ts.interpolated(x,interpolation=interpolation_x, order=order_x,maxbreaks=maxbreaks_x)
  yi <- get.ts.interpolated(y,interpolation=interpolation_y, order=order_y,maxbreaks=maxbreaks_y)
  
  #concert interpolated ts
  zxi <- xi
  zyi <- yi
  #subset interpolated ts
  zxiw <- window(zxi,start = st, end = en)
  zyiw <- window(zyi,start = st, end = en)
  
  ###################################################################
  # calculate magnitude difference between two "real"ts observations
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
  y_mag_forxw <- window(y_mag_forx,start = st, end = en)
  x_mag_foryw <- window(x_mag_fory,start = st, end = en)
  
  #Merge TS
  wc <- merge(zxw,zyw,zxiw,zyiw,x_mag_foryw,y_mag_forxw,zxw[]<-NA,zyw[]<-NA,zyw[]<-NA,zyw[]<-NA)
  names(wc) <- c("x", "y", "xi", "yi","xmag", "ymag","xtimediff", "ytimediff","timediff", "indextimediff")  
  
  #wc[1:200,]
  
  #if only "one-way" only for X-values the y-values have to be related
  if (itype=="one-way"){
    #for all x-observations
    for(i in 1:length(wc$x[!is.na(wc$x)])) {
      #time difference of y value with all x values
      xydiff <- time(wc$x[!is.na(wc$x)][i])-time(wc$y)
      #absolute time difference
      xydiff <- abs(xydiff)  
      #only continue with time difference for values for which y is exisiting
      xydiff[is.na(wc$y)] <- NA
      #check if at least one datapoint is found
      if(length(na.omit(xydiff))>0){
        #wc$magdiff[!is.na(wc$x)][i] <- wc$ymag[!is.na(wc$x)][i]
        wc$xtimediff[!is.na(wc$x)][i] <- wc$x[!is.na(wc$x)][i]
        #index of minimal time difference
        xydiffindex <- which(xydiff == min(na.omit(xydiff)))
        wc$timediff[!is.na(wc$x)][i] <- xydiff[xydiffindex]
        wc$ytimediff[!is.na(wc$x)][i] <- wc$y[xydiffindex]
        wc$indextimediff[!is.na(wc$x)][i] <- as.double(xydiffindex)
      }
    }
    #if "one-way" delete information for all row for which no "real" x-observation exisits
    #wc$ytimediff[is.na(wc$x)] <- NA
    wc$indextimediff[is.na(wc$x)] <- NA
    wc$timediff[is.na(wc$x)] <- NA
    wc$xi[is.na(wc$x)] <- NA
    wc$yi[is.na(wc$x)] <- NA
    wc$xmag[is.na(wc$x)] <- NA
    wc$ymag[is.na(wc$x)] <- NA
  } else if (itype=="two-ways"){
    #for all x-observations
    for(i in 1:length(wc$x[!is.na(wc$x)])) {
      #time difference of y value with all x values
      xydiff <- time(wc$x[!is.na(wc$x)][i])-time(wc$y)
      #absolute time difference
      xydiff <- abs(xydiff)  
      #only continue with time difference for values for which y is exisiting
      xydiff[is.na(wc$y)] <- NA
      #check if at least one datapoint is found
      if(length(na.omit(xydiff))>0){
        #wc$magdiff[!is.na(wc$x)][i] <- wc$ymag[!is.na(wc$x)][i]
        wc$xtimediff[!is.na(wc$x)][i] <- wc$x[!is.na(wc$x)][i]
        #index of minimal time difference
        xydiffindex <- which(xydiff == min(na.omit(xydiff)))
        wc$timediff[!is.na(wc$x)][i] <- xydiff[xydiffindex]
        wc$ytimediff[!is.na(wc$x)][i] <- wc$y[xydiffindex]
        wc$indextimediff[!is.na(wc$x)][i] <- as.double(xydiffindex)
      }
    }
    #for all y-observations
    for(i in 1:length(wc$y[!is.na(wc$y)])) {
      #time(wc[,2][!is.na(wc[,2])][i])-time(wc[,1][!is.na(as.zoo(xw))])
      #time difference of y value with all x values
      xydiff <- time(wc$y[!is.na(wc$y)][i])-time(wc$x)
      #absolute time difference
      xydiff <- abs(xydiff)  
      #only continue with time difference for values for which y is exisiting
      xydiff[is.na(wc$x)] <- NA
      #check if at least one datapoint is found
      if(length(na.omit(xydiff))>0){
        #wc$magdiff[!is.na(wc$y)][i] <- wc$xmag[!is.na(wc$y)][i]
        wc$ytimediff[!is.na(wc$y)][i] <- wc$y[!is.na(wc$y)][i]
        #index of minimal time difference
        xydiffindex <- which(xydiff == min(na.omit(xydiff)))
        wc$timediff[!is.na(wc$y)][i] <- xydiff[xydiffindex]
        wc$xtimediff[!is.na(wc$y)][i] <- wc$x[xydiffindex]
        wc$indextimediff[!is.na(wc$y)][i] <- as.double(xydiffindex)
      }
    }
    #if "one-way" delete information for all row for which no "real" x- or y-observation exisits
    wc$xi[is.na(wc$x)&is.na(wc$y)] <- NA
    wc$yi[is.na(wc$x)&is.na(wc$y)] <- NA
    wc$xmag[is.na(wc$x)&is.na(wc$y)] <- NA
    wc$ymag[is.na(wc$x)&is.na(wc$y)] <- NA
  } else stop("Not a correct interpolation type is selected ('one-way' or 'two-ways') ")
  
  #wc[1:200,]
  
  #if maxdt is provided then delete row information for which timediff > maxdt
  if(!is.null(maxdt)){
    zmaxdt=maxdt/365
    wc$xi[wc$timediff >= zmaxdt] <- NA
    wc$yi[wc$timediff >= zmaxdt] <- NA
    wc$xmag[wc$timediff >= zmaxdt] <- NA
    wc$ymag[wc$timediff >= zmaxdt] <- NA
    wc$xtimediff[wc$timediff >= zmaxdt] <- NA
    wc$ytimediff[wc$timediff >= zmaxdt] <- NA
    #wc$magdiff[wc$timediff >= zmaxdt] <- NA
    #wc$indextimediff[wc$timediff >= zmaxdt] <- NA
    wc$timediff[wc$timediff >= zmaxdt] <- NA
  }
  
  #convert time-difference from zoo to days
  wc[]$timediff <- wc[]$timediff * 365
  return(wc)
}
ts.correlation_A3.plot <- function(x,y,maxdt=NULL,itype=c("one-way","two-ways"),interpolation_x=c("linear","season"),interpolation_y=c("linear","season"),order_x=2,order_y=2,maxbreaks_x=2,maxbreaks_y=2,name_x="NDVI",name_y="HH") {
  
  #Function that returns interpolated TS
  get.ts.interpolated <- function(x,interpolation=c("linear","season"),order=2,maxbreaks=2) {
    if(interpolation=="linear"){
      xi <- na.approx(x)
    } else 
      if (interpolation_x=="season"){
        #get season ts  
        xi <- ts.get_season_ts(x,order=order_x,maxbreaks=maxbreaks_x)
        #substitute points in season ts with real values
        xi[!is.na(as.zoo(x))] <- x[!is.na(as.zoo(x))]
      } else stop("Not a correct interpolation model is selected ('linear' or 'season') ")
    return(xi)
  }
  
  #concert days to zoo
  zx <- as.zoo(x)
  zy <- as.zoo(y)
  #start and end of overlappping period
  st <- max(start(na.omit(zx)),start(na.omit(zy)))
  en <- min(end(na.omit(zx)),end(na.omit(zy)))
  #include one valid observation exceeding the overlapping ts by maxdt days (max and min) 
  if(!is.null(maxdt)){
    zmaxdt=maxdt/365
    st <- st - zmaxdt
    en <- en + zmaxdt
  }
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
  xi <- get.ts.interpolated(x,interpolation=interpolation_x, order=order_x,maxbreaks=maxbreaks_x)
  yi <- get.ts.interpolated(y,interpolation=interpolation_y, order=order_y,maxbreaks=maxbreaks_y)
  
  #concert interpolated ts
  zxi <- xi
  zyi <- yi
  #subset interpolated ts
  zxiw <- window(zxi,start = st, end = en)
  zyiw <- window(zyi,start = st, end = en)
  
  #plot interpolated ts 
  plot(zxiw,xlab="Time",ylab=name_x)
  points(x)
  title(main=paste("Overlapping period [Interpolation: ", interpolation_x,"]",sep=""))  
  plot(zyiw,,xlab="Time",ylab=name_y)
  points(y)
  title(main=paste("Overlapping period [Interpolation: ", interpolation_y,"]",sep=""))  
  plot.2ts(zxiw,zyiw)
  plot.2ts(zxw,zyw)
  plot.2ts(x,y)
  
  ###################################################################
  # calculate magnitude difference between two "real"ts observations
  # all ts points in between two "real" ts observations are assingned
  # with magnitude difference
  glob_st <- min(start(na.omit(zx)),start(na.omit(zy)))
  glob_en <- max(end(na.omit(zx)),end(na.omit(zy)))
  #ye <- as.ts(window(y,start=start(x),end=end(x),extend=TRUE))
  ye <- as.ts(window(y,start=glob_st,end=glob_en,extend=TRUE))
  #get magnitude ts
  y_mag_forx <- abs(get.y_mag_x(x,ye))
  x_mag_fory <- abs(get.y_mag_x(ye,x))
  
  #set magnitude 0 at actual obs
  #y_mag_forx[!is.na(y)] <- 0
  #x_mag_fory[!is.na(x)] <- 0
  #!is.na(x)
  
  #subset magnitude difference ts  
  y_mag_forxw <- window(y_mag_forx,start = st, end = en)
  x_mag_foryw <- window(x_mag_fory,start = st, end = en)
  
  #plot magnitude difference ts
  plot.2ts(zxiw,x_mag_foryw,ylim_x=c(0,1),ylim_y=c(0,1),name_x="NDVI",name_y="Magnitude difference (i,i+1)")
  title(main=paste("X-ts and magnitude difference between obs (used for y-ts obs)"),cex.main=0.9)
  plot.2ts(zyiw,y_mag_forxw,name_x="HHHV",name_y="Magnitude difference (i,i+1)")
  title(main=paste("Y-ts and magnitude difference between obs (used for x-ts obs)"),cex.main=0.9)
  
  #Merge TS
  wc <- merge(zxw,zyw,zxiw,zyiw,x_mag_foryw,y_mag_forxw,zxw[]<-NA,zyw[]<-NA,zyw[]<-NA,zyw[]<-NA)
  names(wc) <- c("x", "y", "xi", "yi","xmag", "ymag","xtimediff", "ytimediff","timediff", "indextimediff")  
  head(wc) 
  wc[1:300,]
  
  #if only "one-way" only for X-values the y-values have to be related
  if (itype=="one-way"){
    #for all x-observations
    for(i in 1:length(wc$x[!is.na(wc$x)])) {
      #time difference of y value with all x values
      xydiff <- time(wc$x[!is.na(wc$x)][i])-time(wc$y)
      #absolute time difference
      xydiff <- abs(xydiff)  
      #only continue with time difference for values for which y is exisiting
      xydiff[is.na(wc$y)] <- NA
      #check if at least one datapoint is found
      if(length(na.omit(xydiff))>0){
        #wc$magdiff[!is.na(wc$x)][i] <- wc$ymag[!is.na(wc$x)][i]
        wc$xtimediff[!is.na(wc$x)][i] <- wc$x[!is.na(wc$x)][i]
        #index of minimal time difference
        xydiffindex <- which(xydiff == min(na.omit(xydiff)))
        wc$timediff[!is.na(wc$x)][i] <- xydiff[xydiffindex]
        wc$ytimediff[!is.na(wc$x)][i] <- wc$y[xydiffindex]
        wc$indextimediff[!is.na(wc$x)][i] <- as.double(xydiffindex)
      }
    }
    #if "one-way" delete information for all row for which no "real" x-observation exisits
    #wc$ytimediff[is.na(wc$x)] <- NA
    wc$indextimediff[is.na(wc$x)] <- NA
    wc$timediff[is.na(wc$x)] <- NA
    wc$xi[is.na(wc$x)] <- NA
    wc$yi[is.na(wc$x)] <- NA
    wc$xmag[is.na(wc$x)] <- NA
    wc$ymag[is.na(wc$x)] <- NA
  } else if (itype=="two-ways"){
    #for all x-observations
    for(i in 1:length(wc$x[!is.na(wc$x)])) {
      #time difference of y value with all x values
      xydiff <- time(wc$x[!is.na(wc$x)][i])-time(wc$y)
      #absolute time difference
      xydiff <- abs(xydiff)  
      #only continue with time difference for values for which y is exisiting
      xydiff[is.na(wc$y)] <- NA
      #check if at least one datapoint is found
      if(length(na.omit(xydiff))>0){
        #wc$magdiff[!is.na(wc$x)][i] <- wc$ymag[!is.na(wc$x)][i]
        wc$xtimediff[!is.na(wc$x)][i] <- wc$x[!is.na(wc$x)][i]
        #index of minimal time difference
        xydiffindex <- which(xydiff == min(na.omit(xydiff)))
        wc$timediff[!is.na(wc$x)][i] <- xydiff[xydiffindex]
        wc$ytimediff[!is.na(wc$x)][i] <- wc$y[xydiffindex]
        wc$indextimediff[!is.na(wc$x)][i] <- as.double(xydiffindex)
      }
    }
    #for all y-observations
    for(i in 1:length(wc$y[!is.na(wc$y)])) {
      #time(wc[,2][!is.na(wc[,2])][i])-time(wc[,1][!is.na(as.zoo(xw))])
      #time difference of y value with all x values
      xydiff <- time(wc$y[!is.na(wc$y)][i])-time(wc$x)
      #absolute time difference
      xydiff <- abs(xydiff)  
      #only continue with time difference for values for which y is exisiting
      xydiff[is.na(wc$x)] <- NA
      #check if at least one datapoint is found
      if(length(na.omit(xydiff))>0){
        #wc$magdiff[!is.na(wc$y)][i] <- wc$xmag[!is.na(wc$y)][i]
        wc$ytimediff[!is.na(wc$y)][i] <- wc$y[!is.na(wc$y)][i]
        #index of minimal time difference
        xydiffindex <- which(xydiff == min(na.omit(xydiff)))
        wc$timediff[!is.na(wc$y)][i] <- xydiff[xydiffindex]
        wc$xtimediff[!is.na(wc$y)][i] <- wc$x[xydiffindex]
        wc$indextimediff[!is.na(wc$y)][i] <- as.double(xydiffindex)
      }
    }
    #if "one-way" delete information for all row for which no "real" x- or y-observation exisits
    wc$xi[is.na(wc$x)&is.na(wc$y)] <- NA
    wc$yi[is.na(wc$x)&is.na(wc$y)] <- NA
    wc$xmag[is.na(wc$x)&is.na(wc$y)] <- NA
    wc$ymag[is.na(wc$x)&is.na(wc$y)] <- NA
  } else stop("Not a correct interpolation type is selected ('one-way' or 'two-ways') ")
  
  wc[1:200,]
  plot.2ts(wc$x,wc$y)
  
  
  #if maxdt is provided then delete row information for which timediff > maxdt
  if(!is.null(maxdt)){
    zmaxdt=maxdt/365
    wc$xi[wc$timediff >= zmaxdt] <- NA
    wc$yi[wc$timediff >= zmaxdt] <- NA
    wc$xmag[wc$timediff >= zmaxdt] <- NA
    wc$ymag[wc$timediff >= zmaxdt] <- NA
    wc$xtimediff[wc$timediff >= zmaxdt] <- NA
    wc$ytimediff[wc$timediff >= zmaxdt] <- NA
    #wc$magdiff[wc$timediff >= zmaxdt] <- NA
    #wc$indextimediff[wc$timediff >= zmaxdt] <- NA
    wc$timediff[wc$timediff >= zmaxdt] <- NA
  }
  
  #plots
  plot.2ts(na.omit(wc$xi),na.omit(wc$yi),name_x="Interpolated (x)",name_y="Interpolated (y)")
  calc.Rsquare.plot(na.omit(wc$xi),na.omit(wc$yi),name_x="Interpolated (x)",name_y="Interpolated (y)")
  #plot.2ts(na.omit(wc$xtimediff)),na.omit(wc$ytimediff)),"maxdt (x)",name_y="maxdt (y)")
  #calc.Rsquare.plot(na.omit(wc$xtimediff),na.omit(wc$ytimediff),name_x="maxdt (x)",name_y="maxdt (y)")
  
  #convert time-difference from zoo to days
  wc[]$timediff <- wc[]$timediff * 365
  #Print summary
  print(wc[!is.na(wc$x)])
  print(wc[!is.na(wc$y)])
  # overview
  return(wc)
}

#!! Layover and shadow area need to be masked first!!
calc.brick_ts.correlation_A3_opt_wt <- function(rx,ry,dx,dy,maxdt=NULL,itype=itype,interpolation_x=c('linear','season'),interpolation_y=c('linear','season'),order_x=1,order_y=1,maxbreaks_x = 2, maxbreaks_y = 2,max_wt_exp=10,steps=1,cores=1){
  #Output:
  #Layer1: Qualit flag
  ##QF = 1 when not: (i) ts overlap and at least 2 obs are overlapping in each ts
  ##QF = 2 if ts returns less than 2 correlation pairs
  #Layer2 (rg): range of overlapping time period of xts over yts in years
  calc_c <- function(v){
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[length(dx)+1:length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    
    #Quality flag
    #QF == 1 when not: (i) ts overlap and at least 2 obs are overlapping in each ts
    #QF == 2 if ts returns less than 2 correlation pairs
    QF <- 0
    #Check if (i) ts overlap and at least 2 points are overlapping in each ts, (ii) ts have at least 2 points
    #if(isTRUE(ts.overlap(ddx,ddy-ddz)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy-ddz)))>1){
    if(min(ts.overlap_nr(ddx,ddy)[]) > 1){
      wc <- ts.correlation_A3(ddx,ddy,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      xi <- wc$xi
      yi <- wc$yi
      #linear regression function
      #if pairs < 2 regression does not work
      if(length(na.omit(as.zoo(xi))) >= 2){
        wmag <- get.magweight(wc$x,wc$y,wc$xmag,wc$ymag,itype="two-ways")
        a <- opt_wt(xi,yi,wmag,max_wt_exp=max_wt_exp,steps=steps,order=1)
        opt_pVal <- as.double(a$pVal)
        wt_exp <- a$wt_exp
        
        vxi <- as.vector(xi)
        vyi <- as.vector(yi)
        fit=lm(na.omit(vyi)~poly(na.omit(vxi),1),weight=wmag^wt_exp)
        pVal <- as.double(anova(fit)$'Pr(>F)'[1])
        rsquare <- summary(fit)$r.squared
        rsquare_lci <- as.double(CIr(r=rsquare, n = length(na.omit(as.zoo(xi))), level = .90))[1]
        rg <- range(na.omit(index(na.omit(wc$x))))[2]-range(na.omit(index(na.omit(wc$x))))[1]
        n = length(na.omit(as.zoo(xi)))
      } else {
        rsquare <- NA
        pVal <- NA
        opt_pVal <- NA
        rsquare_lci <- NA
        n <-NA
        wt_exp <- NA
        rg <- NA
        #QF == 2 if ts returns less than 2 correlation pairs
        QF <- 2
      }
      #rsquare <- as.double(rsquare)
    } else {
      rsquare <- NA
      pVal <- NA
      opt_pVal <- NA
      #rsquare_wtimediff <- NA
      # rsquare_wcom <- NA
      rsquare_lci <- NA
      n <-NA
      wt_exp <- NA
      rg <- NA
      #QF == 2 if ts returns less than 2 correlation pairs
      QF <- 1
    }
    #combine vectors
    c <- c(QF,rg,n,rsquare,rsquare_lci,wt_exp,pVal,opt_pVal)
    #c <- c(QF,rg,rsquare,rsquare_lci)
    return(c)
  }
  
  fun <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- matrix(NA, length(i), 8)
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1, calc_c))
    }
    return(res)
  }
  
  #combine input raster bricks
  b <- brick(addLayer(rx,ry))
  
  if(cores==1){
    print("start")
    x <- mc.calc(b,fun=fun)
  }else {
    print("start")
    x <- mc.calc(b, fun=fun,mc.cores=cores)
  }
  names(x) <- c("quality flag","overlap years","n","r-square","LCIr-square90","wt_exp","p-value","p-value (opt)")
  return(x)
}
calc.brick_ts.correlation_A3 <- function(rx,ry,dx,dy,wt_exp=0,maxdt=NULL,itype=itype,interpolation_x=c('linear','season'),interpolation_y=c('linear','season'),order_x=1,order_y=1,maxbreaks_x = 2, maxbreaks_y = 2,cores=1){
  #Output:
  #Layer1: Qualit flag
  ##QF = 1 when not: (i) ts overlap and at least 2 obs are overlapping in each ts
  ##QF = 2 if ts returns less than 2 correlation pairs
  calc_c_hhhv <- function(v){
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[length(dx)+1:length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    
    #Quality flag
    #QF == 1 when not: (i) ts overlap and at least 2 obs are overlapping in each ts
    #QF == 2 if ts returns less than 2 correlation pairs
    QF <- 0
    #Check if (i) ts overlap and at least 2 points are overlapping in each ts, (ii) ts have at least 2 points
    #if(isTRUE(ts.overlap(ddx,ddy-ddz)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy-ddz)))>1){
    if(min(ts.overlap_nr(ddx,ddy)[]) > 1){
      wc <- ts.correlation_A3(ddx,ddy,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      xi <- wc$xi
      yi <- wc$yi
      #linear regression function
      #if pairs < 2 regression does not work
      if(length(na.omit(as.zoo(xi))) >= 2){
        wmag <- get.magweight(wc$x,wc$y,wc$xmag,wc$ymag,itype="two-ways")
        #wcom <- wmag + wtimediff
        vxi <- as.vector(xi)
        vyi <- as.vector(yi)
        fit=lm(na.omit(vyi)~poly(na.omit(vxi),1),weight=wmag^wt_exp)
        pVal <- as.double(anova(fit)$'Pr(>F)'[1])
        rsquare <- summary(fit)$r.squared
        rsquare_lci <- as.double(CIr(r=rsquare, n = length(na.omit(as.zoo(xi))), level = .90))[1]
        rg <- range(na.omit(index(na.omit(wc$x))))[2]-range(na.omit(index(na.omit(wc$x))))[1]
        n <- length(na.omit(as.zoo(xi)))
      } else {
        rsquare <- NA
        pVal <- NA
        rsquare_lci <- NA
        n <-NA
        rg <- NA
        #QF == 2 if ts returns less than 2 correlation pairs
        QF <- 2
      }
      #rsquare <- as.double(rsquare)
    } else {
      rsquare <- NA
      pVal <- NA
      rsquare_lci <- NA
      n <-NA
      rg <- NA
      #QF = 1 when not: (i) ts overlap, (ii) ts have at least 2 points
      QF <- 1
    }
    #combine vectors
    c <- c(QF,rg,n,rsquare,rsquare_lci,wt_exp,pVal)
    return(c)
  }
  
  fun <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- matrix(NA, length(i), 7)
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1, calc_c_hhhv))
    }
    return(res)
  }
  
  #combine input raster bricks
  b <- brick(addLayer(rx,ry))
  
  if(cores==1){
    print("start")
    x <- mc.calc(b,fun=fun)
  }else {
    print("start")
    x <- mc.calc(b, fun=fun,mc.cores=cores)
  }
  names(x) <- c("quality flag","overlap years","n","r-square","LCIr-square90","wt_exp","p-value")
  return(x)
}

#get brick for entire sa
get.brick_ts.correlation_A3 <- function(rx,ry,dx,dy,maxdt=NULL,itype=itype,interpolation_x=c('linear','season'),interpolation_y=c('linear','season'),order_x=1,order_y=1,maxbreaks_x = 2, maxbreaks_y = 2,cores=1){
  #calc length overlap
  calc.length_overlap <- function(){
    ddx <- bfastts(as.vector(rx[1]),dx,type=c("irregular"))
    ddy <- bfastts(as.vector(ry[1]),dy,type=c("irregular"))
    #ddz <- bfastts(as.vector(rz[1]),dz,type=c("irregular"))
    #convert to zoo
    zx <- as.zoo(ddx)
    zy <- as.zoo(ddy)
    st <- max(start(zx),start(zy))
    en <- min(end(zx),end(zy))
    #extend overlapping period
    #include one valid observation exceeding the overlapping ts by maxdt days (max and min) 
    #merge ts based on the longest ts to extend the shorter one to the length of the longer
    #necessary since the command "window" does not exceed ts and this would result in two differently long ts
    #subset ts and merge
    zxw <- window(zy,start = st, end = en)
    length_overlap <- length(zxw)
    return(length_overlap)
  }
  
  #get x ts for brick cell
  get_x <- function(v){
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[length(dx)+1:length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    #if(isTRUE(ts.overlap(ddx,ddy)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy)))>1){
    if(min(ts.overlap_nr(ddx,ddy)[]) > 1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cbx <- wc$x
      #cb <- ts.correlation(ddx,ddy,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      #cbx <- cb[,1]
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddx_new <- rep(NA, calc.length_overlap())
      ddx_new[1:length(cbx)] <- cbx[]
    } else {
      ddx_new <- rep(NA, calc.length_overlap())
    }
    return(ddx_new)
  }
  
  #get y ts for brick cell
  get_y <- function(v){
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[length(dx)+1:length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    if(isTRUE(ts.overlap(ddx,ddy)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy)))>1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cby <- wc$y
      #cb <- ts.correlation(ddx,ddy,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y) 
      #cby <- cb[,2]
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddy_new <- rep(NA, calc.length_overlap())
      ddy_new[1:length(cby)] <- cby[]
    } else {
      ddy_new <- rep(NA, calc.length_overlap())
    }
    return(ddy_new)
  }
  
  #get x ts for brick cell
  get_xi <- function(v){
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[length(dx)+1:length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    if(isTRUE(ts.overlap(ddx,ddy)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy)))>1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cbx <- wc$xi
      #cb <- ts.correlation(ddx,ddy,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      #cbx <- cb[,1]
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddx_new <- rep(NA, calc.length_overlap())
      ddx_new[1:length(cbx)] <- cbx[]
    } else {
      ddx_new <- rep(NA, calc.length_overlap())
    }
    return(ddx_new)
  }
  
  #get y ts for brick cell
  get_yi <- function(v){
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[length(dx)+1:length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    if(isTRUE(ts.overlap(ddx,ddy)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy)))>1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cby <- wc$yi
      #cb <- ts.correlation(ddx,ddy,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y) 
      #cby <- cb[,2]
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddy_new <- rep(NA, calc.length_overlap())
      ddy_new[1:length(cby)] <- cby[]
    } else {
      ddy_new <- rep(NA, calc.length_overlap())
    }
    return(ddy_new)
  }
  
  #get x ts for brick cell
  get_xmag <- function(v){
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[length(dx)+1:length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    if(isTRUE(ts.overlap(ddx,ddy)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy)))>1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cbx <- wc$xmag
      #cb <- ts.correlation(ddx,ddy,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      #cbx <- cb[,1]
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddx_new <- rep(NA, calc.length_overlap())
      ddx_new[1:length(cbx)] <- cbx[]
    } else {
      ddx_new <- rep(NA, calc.length_overlap())
    }
    return(ddx_new)
  }
  
  #get x ts for brick cell
  get_ymag <- function(v){
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[length(dx)+1:length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    if(isTRUE(ts.overlap(ddx,ddy)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy)))>1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cbx <- wc$ymag
      #cb <- ts.correlation(ddx,ddy,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      #cbx <- cb[,1]
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddx_new <- rep(NA, calc.length_overlap())
      ddx_new[1:length(cbx)] <- cbx[]
    } else {
      ddx_new <- rep(NA, calc.length_overlap())
    }
    return(ddx_new)
  }
  
  
  funx <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_x))        
    }
    return(res)
  }
  
  funy <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_y))        
    }
    return(res)
  }
  
  funxi <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_xi))        
    }
    return(res)
  }
  
  funyi <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_yi))        
    }
    return(res)
  }
  
  fun_xmag <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_xmag))        
    }
    return(res)
  }
  
  fun_ymag <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_ymag))        
    }
    return(res)
  }
  
  
  
  b <- brick(addLayer(rx,ry))
  
  if(cores==1){
    x <- calc(b, fun=funx)
    y <- calc(b, fun=funy)
    xi <- calc(b, fun=funxi)
    yi <- calc(b, fun=funyi)
    xmag <- calc(b, fun=fun_xmag)
    ymag <- calc(b, fun=fun_ymag)
  }else {
    print(sprintf("start"))
    x <- mc.calc(b, fun=funx,mc.cores=cores)
    print(sprintf("past x"))
    y <- mc.calc(b, fun=funy,mc.cores=cores)
    print(sprintf("past y"))
    xi <- mc.calc(b, fun=funxi,mc.cores=cores)
    print(sprintf("past xi"))
    yi <- mc.calc(b, fun=funyi,mc.cores=cores)
    print(sprintf("past yi"))
    xmag <- mc.calc(b, fun=fun_xmag,mc.cores=cores)
    print(sprintf("past xmag"))
    ymag <- mc.calc(b, fun=fun_ymag,mc.cores=cores)
    print(sprintf("end"))
  }
  return(cbind(x,y,xi,yi,xmag,ymag))
}
get.brick_ts.correlation_A3_hhhv <- function(rx,ry,rz,dx,dy,dz,maxdt=NULL,itype=itype,interpolation_x=c('linear','season'),interpolation_y=c('linear','season'),order_x=1,order_y=1,maxbreaks_x = 2, maxbreaks_y = 2,cores=1){
  #calculate length of overlapping period
  #overlapping period is extended by +- maxdt to include possible correlation partners
  #in the range of maxdt
  calc.length_overlap_hhhv <- function(){
    ddx <- bfastts(as.vector(rx[1]),dx,type=c("irregular"))
    ddy <- bfastts(as.vector(ry[1]),dy,type=c("irregular"))
    ddz <- bfastts(as.vector(rz[1]),dz,type=c("irregular"))
    #convert to zoo
    zx <- as.zoo(ddx)
    zy <- as.zoo(ddy-ddz)
    st <- max(start(zx),start(zy))
    en <- min(end(zx),end(zy))
    #extend overlapping period
    #include one valid observation exceeding the overlapping ts by maxdt days (max and min) 
    #merge ts based on the longest ts to extend the shorter one to the length of the longer
    #necessary since the command "window" does not exceed ts and this would result in two differently long ts
    #subset ts and merge
    zxw <- window(zy,start = st, end = en)
    length_overlap <- length(zxw)
    return(length_overlap)
  }
  
  #get x ts for brick cell
  get_x_hhhv <- function(v){
    #split up combined input ts
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[(length(dx)+1):(length(dx)+length(dy))]
    z <- v[(length(dx)+length(dy)+1):length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    ddz <- bfastts(z,dz,type=c("irregular"))
    #get correlation pairs 
    #cb <- ts.correlation(ddx,ddy-ddz,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
    #return(cb[,1])
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    #if(isTRUE(ts.overlap(ddx,ddy-ddz)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy-ddz)))>1){
    if(min(ts.overlap_nr(ddx,ddy-ddz)[]) > 1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy-ddz,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cbx <- wc$x
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddx_new <- rep(NA, calc.length_overlap_hhhv())
      ddx_new[1:length(cbx)] <- cbx[]
    } else {
      ddx_new <- rep(NA, calc.length_overlap_hhhv())
    }
    return(ddx_new)
    #return(sum(na.omit(y)))
  }
  #get y ts for brick cell
  get_y_hhhv <- function(v){
    x <- v[1:length(dx)]
    y <- v[(length(dx)+1):(length(dx)+length(dy))]
    z <- v[(length(dx)+length(dy)+1):length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    ddz <- bfastts(z,dz,type=c("irregular"))
    #get correlation pairs 
    #cb <- ts.correlation(ddx,ddy-ddz,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
    #return(cb[,2])
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    #if(isTRUE(ts.overlap(ddx,ddy-ddz)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy-ddz)))>1){
    if(min(ts.overlap_nr(ddx,ddy-ddz)[]) > 1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy-ddz,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cby <- wc$y
      #cb <- ts.correlation(ddx,ddy-ddz,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y) 
      #cby <- cb[,2]
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddy_new <- rep(NA, calc.length_overlap_hhhv())
      ddy_new[1:length(cby)] <- cby[]
    } else {
      ddy_new <- rep(NA, calc.length_overlap_hhhv())
    }
    return(ddy_new)
  }
  #get x ts for brick cell
  get_xi_hhhv <- function(v){
    #split up combined input ts
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[(length(dx)+1):(length(dx)+length(dy))]
    z <- v[(length(dx)+length(dy)+1):length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    ddz <- bfastts(z,dz,type=c("irregular"))
    #get correlation pairs 
    #cb <- ts.correlation(ddx,ddy-ddz,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
    #return(cb[,1])
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    #if(isTRUE(ts.overlap(ddx,ddy-ddz)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy-ddz)))>1){
    if(min(ts.overlap_nr(ddx,ddy-ddz)[]) > 1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy-ddz,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cbx <- wc$xi
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddx_new <- rep(NA, calc.length_overlap_hhhv())
      ddx_new[1:length(cbx)] <- cbx[]
    } else {
      ddx_new <- rep(NA, calc.length_overlap_hhhv())
    }
    return(ddx_new)
    #return(sum(na.omit(y)))
  }
  #get y ts for brick cell
  get_yi_hhhv <- function(v){
    x <- v[1:length(dx)]
    y <- v[(length(dx)+1):(length(dx)+length(dy))]
    z <- v[(length(dx)+length(dy)+1):length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    ddz <- bfastts(z,dz,type=c("irregular"))
    #get correlation pairs 
    #cb <- ts.correlation(ddx,ddy-ddz,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
    #return(cb[,2])
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    #if(isTRUE(ts.overlap(ddx,ddy-ddz)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy-ddz)))>1){
    if(min(ts.overlap_nr(ddx,ddy-ddz)[]) > 1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy-ddz,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cby <- wc$yi
      #cb <- ts.correlation(ddx,ddy-ddz,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y) 
      #cby <- cb[,2]
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddy_new <- rep(NA, calc.length_overlap_hhhv())
      ddy_new[1:length(cby)] <- cby[]
    } else {
      ddy_new <- rep(NA, calc.length_overlap_hhhv())
    }
    return(ddy_new)
  }
  
  #get y ts for brick cell
  get_xmag_hhhv <- function(v){
    x <- v[1:length(dx)]
    y <- v[(length(dx)+1):(length(dx)+length(dy))]
    z <- v[(length(dx)+length(dy)+1):length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    ddz <- bfastts(z,dz,type=c("irregular"))
    #get correlation pairs 
    #cb <- ts.correlation(ddx,ddy-ddz,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
    #return(cb[,2])
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    #if(isTRUE(ts.overlap(ddx,ddy-ddz)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy-ddz)))>1){
    if(min(ts.overlap_nr(ddx,ddy-ddz)[]) > 1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy-ddz,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cby <- wc$xmag
      #cb <- ts.correlation(ddx,ddy-ddz,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y) 
      #cby <- cb[,2]
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddy_new <- rep(NA, calc.length_overlap_hhhv())
      ddy_new[1:length(cby)] <- cby[]
    } else {
      ddy_new <- rep(NA, calc.length_overlap_hhhv())
    }
    return(ddy_new)
  }
  #get y ts for brick cell
  get_ymag_hhhv <- function(v){
    x <- v[1:length(dx)]
    y <- v[(length(dx)+1):(length(dx)+length(dy))]
    z <- v[(length(dx)+length(dy)+1):length(v)]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    ddz <- bfastts(z,dz,type=c("irregular"))
    #get correlation pairs 
    #cb <- ts.correlation(ddx,ddy-ddz,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
    #return(cb[,2])
    #Check if (i) ts overlap, (ii) ts have at least 2 points
    #if(isTRUE(ts.overlap(ddx,ddy-ddz)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy-ddz)))>1){
    if(min(ts.overlap_nr(ddx,ddy-ddz)[]) > 1){
      #get correlation pairs 
      wc <- ts.correlation_A3(ddx,ddy-ddz,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      cby <- wc$ymag
      #cb <- ts.correlation(ddx,ddy-ddz,itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y) 
      #cby <- cb[,2]
      #the output length of ts.correlation is varying, but we need a fixed length (length_overlap)
      #thus, (i) create NA vector with length overlap, (ii) replace the output
      ddy_new <- rep(NA, calc.length_overlap_hhhv())
      ddy_new[1:length(cby)] <- cby[]
    } else {
      ddy_new <- rep(NA, calc.length_overlap_hhhv())
    }
    return(ddy_new)
  }
  
  funx <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap_hhhv()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_x_hhhv))        
    }
    return(res)
  }
  funy <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap_hhhv()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_y_hhhv))        
    }
    return(res)
  }
  funxi <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap_hhhv()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_xi_hhhv))        
    }
    return(res)
  }
  funyi <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap_hhhv()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_yi_hhhv))        
    }
    return(res)
  }
  fun_xmag <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap_hhhv()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_xmag_hhhv))        
    }
    return(res)
  }
  fun_ymag <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- array(NA,c(length(i),calc.length_overlap_hhhv()))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      res[i, ] <- t(apply(y[i, ], 1,get_ymag_hhhv))        
    }
    return(res)
  }
  
  b <- brick(addLayer(rx,ry,rz))
  if(cores==1){
    x <- calc(b, fun=funx)
    y <- calc(b, fun=funy)
    xi <- calc(b, fun=funxi)
    yi <- calc(b, fun=funyi)
    xmag <- calc(b, fun=fun_xmag)
    ymag <- calc(b, fun=fun_ymag)
    td <- calc(b, fun=fun_td)
  }else {
    print(sprintf("start"))
    x <- mc.calc(b, fun=funx,mc.cores=cores)
    print(sprintf("past x"))
    y <- mc.calc(b, fun=funy,mc.cores=cores)
    print(sprintf("past y"))
    xi <- mc.calc(b, fun=funxi,mc.cores=cores)
    print(sprintf("past xi"))
    yi <- mc.calc(b, fun=funyi,mc.cores=cores)
    print(sprintf("past yi"))
    xmag <- mc.calc(b, fun=fun_xmag,mc.cores=cores)
    print(sprintf("past xmag"))
    ymag <- mc.calc(b, fun=fun_ymag,mc.cores=cores)
    print(sprintf("end"))
  }
  return(cbind(x,y,xi,yi,xmag,ymag))
}

##Fusion
ts.fusion.correlation_A3 <- function(ddx,ddy,wt_type = c("none", "magnitude"),wt_exp=0,wt_opt=FALSE,wt_opt_steps=1,wt_opt_max_exp=5,maxdt = NULL,itype=c("one-way","two-ways"),interpolation_x=c('linear','season'),interpolation_y=c('linear','season'),order_x=1,order_y=1,maxbreaks_x=5,maxbreaks_y=5){
  # itype: one-way:  Interpolation for all data point for which Xexisits
  # itype: tow-ways:  Interpolation for all data point for which X or Y exisits
  wc <- ts.correlation_A3(ddx,ddy,itype=itype,maxdt=maxdt,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_x)
  xi <- as.zoo(wc[,3])
  yi <- as.zoo(wc[,4])
  
  if(length(na.omit(xi)) < 2) {
    ddx[]<-NA
  } else {
    if(wt_type == "none"){
      #no weight
      #Regression: 
      #lm(Y ~ model); where Y is the dependent variable to be predicted by model "poly(x,1)"
      #y <- x and x <- y
      y <- as.vector(na.omit(xi))
      x <- as.vector(na.omit(yi))
      #lm model
    } else {
      #calculate weights
      if(wt_type =="magnitude"){wt <- get.magweight(wc$x,wc$y,wc$xmag,wc$ymag,itype=itype)}
      if(wt_opt == TRUE){
        a <- opt_wt(xi,yi,wt,max_wt_exp=wt_opt_max_exp,steps=wt_opt_steps,order=1)
        wt_exp <- a$wt_exp
      }
      
      #Regression: 
      #lm(Y ~ model); where Y is the dependent variable to be predicted by model "poly(x,1)"
      #y <- x and x <- y
      y <- as.vector(na.omit(xi))
      x <- as.vector(na.omit(yi))
      #lm model
      fit=lm(y~poly(x,1),weight=wt^wt_exp)
    }
    
    #prediction
    new_data <- na.omit(as.zoo(ddy))
    pre <- predict(fit, data.frame(x=new_data))
    
    #replace predicted x values in x ts at position of y data that was used as input for prediction
    ddx[match(index(as.zoo(na.omit(as.zoo(ddy)))),index(as.zoo(ddx)))] <- as.zoo(pre)
    
  }
  
  return(ddx)
}
ts.fusion.correlation_A3.plot <- function(ddx,ddy,wt_type = c("none", "magnitude"),wt_exp=0,wt_opt=FALSE,wt_opt_steps=1,wt_opt_max_exp=5,maxdt = NULL,itype=c("one-way","two-ways"),interpolation_x=c('linear','season'),interpolation_y=c('linear','season'),order_x=1,order_y=1,maxbreaks_x=5,maxbreaks_y=5,name_x="",name_y=""){
  # itype: one-way:  Interpolation for all data point for which Xexisits
  # itype: tow-ways:  Interpolation for all data point for which X or Y exisits
  wc <- ts.correlation_A3(ddx,ddy,itype=itype,maxdt=maxdt,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_x)
  xi <- as.zoo(wc[,3])
  yi <- as.zoo(wc[,4])
  
  if(length(na.omit(xi)) < 2) {
    ddx[]<-NA
  } else {
    #plot
    plot.2ts(ddx,ddy,name_x=name_x,name_y=name_y,points=TRUE, title="Original TS")
    plot.2ts(xi,yi,name_x=name_x,name_y=name_y,points=TRUE, title="Origianal X TS related Y TS")
    
    if(wt_type == "none"){
      #no weight
      #Regression: 
      #lm(Y ~ model); where Y is the dependent variable to be predicted by model "poly(x,1)"
      #y <- x and x <- y
      y <- as.vector(na.omit(xi))
      x <- as.vector(na.omit(yi))
      #lm model
    } else {
      #calculate weights
      if(wt_type =="magnitude"){wt <- get.magweight(wc$x,wc$y,wc$xmag,wc$ymag,itype=itype)}
      if(wt_opt == TRUE){
        a <- opt_wt.plot(xi,yi,wt,max_wt_exp=wt_opt_max_exp,steps=wt_opt_steps,order=1)
        wt_exp <- a$wt_exp
        print("OPT")
      }
      
      print(paste("wt_exp: ",wt_exp))
      calc.Rsquare_weight.plot(xi,yi,weight=wt^wt_exp,order=1,name_x,name_y)
      #Regression: 
      #lm(Y ~ model); where Y is the dependent variable to be predicted by model "poly(x,1)"
      #y <- x and x <- y
      y <- as.vector(na.omit(xi))
      x <- as.vector(na.omit(yi))
      #lm model
      fit=lm(y~poly(x,1),weight=wt^wt_exp)
    }
    
    #prediction
    new_data <- na.omit(as.zoo(ddy))
    pre <- predict(fit, data.frame(x=new_data))
    
    #replace predicted x values in x ts at position of y data that was used as input for prediction
    ddx[match(index(as.zoo(na.omit(as.zoo(ddy)))),index(as.zoo(ddx)))] <- as.zoo(pre)
    
  }
  
  return(ddx)
}
get.brick_ts.fusion.correlation_A3 <- function(rx,ry,dx,dy,itype=c("one-way","two-ways"),wt_type = c("none", "magnitude"),wt_exp=0,wt_opt=FALSE,wt_opt_steps=1,wt_opt_max_exp=5,maxdt = NULL,interpolation_x=c('linear','season'),interpolation_y=c('linear','season'),order_x=1,order_y=1,maxbreaks_x = 2, maxbreaks_y = 2,cores=1){
  print(paste("itype: ",itype))
  print(paste("wt_type: ",wt_type))
  print(paste("wt_exp: ",wt_exp))
  print(paste("wt_opt: ",wt_opt))
  print(paste("wt_opt_steps: ",wt_opt_steps))
  print(paste("wt_opt_max_exp: ",wt_opt_max_exp))
  #get x ts for brick cell
  get_x <- function(v){
    #split up combined input ts
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[(length(dx)+1):(length(dx)+length(v))]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    #get correlation pairs 
    #cb <- get.ts.correlation(ddx,ddy-ddz,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
    #Check if (i) TS overlap (ii) TS at least 2 valid data points
    #if(isTRUE(ts.overlap(ddx,ddy)) && length(na.omit(as.zoo(ddx)))>1 && length(na.omit(as.zoo(ddy)))>1){
    if(min(ts.overlap_nr(ddx,ddy)[]) > 1){
      ddx_new <- ts.fusion.correlation_A3(ddx,ddy,itype=itype,wt_type=wt_type,wt_exp=wt_exp,wt_opt=wt_opt,wt_opt_steps=wt_opt_steps,wt_opt_max_exp=wt_opt_max_exp,maxdt=maxdt,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      #ddx_new <- ts.fusion.correlation_A3(ddx,ddy,itype=itype,wt_type=wt_type,wt_exp=wt_exp,maxdt=maxdt,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      
    } else {
      ddx_new <- ddx
      ddx_new[] <- NA
    }
    return(ddx_new)
    #return(sum(na.omit(y)))
  }
  
  
  funx <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    length_overlap <- length(bfastts(as.vector(rx[1]),dx,type=c("irregular")))
    res <- array(NA,c(length(i),length_overlap))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      # res[i, ] <- t(apply(y[i, ], l, calc_c))
      res[i, ] <- t(apply(y[i, ], 1,get_x))        
    }
    return(res)
  }
  
  b <- brick(addLayer(rx,ry))
  if(cores==1){
    x <- calc(b, fun=funx)
  }else {
    x <- mc.calc(b, fun=funx, mc.cores=cores)
  }
  #names(x) <- index(bfastts(as.vector(rx[1]),dx,type=c("irregular")))
  return(x)
}
get.brick_ts.fusion.correlation_flags_A3 <- function(rx,ry,dx,dy,itype=c("one-way","two-ways"),wt_type = c("none", "magnitude"),wt_exp=0,wt_opt=FALSE,wt_opt_steps=1,wt_opt_max_exp=5,maxdt = NULL,interpolation_x=c('linear','season'),interpolation_y=c('linear','season'),order_x=1,order_y=1,maxbreaks_x = 2, maxbreaks_y = 2,cores=1){
  #getfusion brick with flags, wt_exp gets calculated and comitted to ts.fusion.correlation_A3
  #therfore wt_opt is set to always FALSE!
  print(paste("itype: ",itype))
  print(paste("wt_type: ",wt_type))
  print(paste("wt_exp: ",wt_exp))
  print(paste("wt_opt: ",wt_opt))
  print(paste("wt_opt_steps: ",wt_opt_steps))
  print(paste("wt_opt_max_exp: ",wt_opt_max_exp))
  #get x ts for brick cell
  get_x <- function(v){
    #split up combined input ts
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[(length(dx)+1):(length(dx)+length(v))]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    #get correlation pairs 
    #cb <- get.ts.correlation(ddx,ddy-ddz,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
    #Check if (i) TS overlap (ii) TS at least 2 valid data points
    
    QF <- 0    
    if(min(ts.overlap_nr(ddx,ddy)[]) > 1){
      wc <- ts.correlation_A3(ddx,ddy,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      xi <- wc$xi
      yi <- wc$yi
      #ddx_new <- ts.fusion.correlation_A3(ddx,ddy,itype=itype,wt_type=wt_type,wt_exp=wt_exp,maxdt=maxdt,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      if(length(na.omit(as.zoo(xi))) >= 2){
        wmag <- get.magweight(wc$x,wc$y,wc$xmag,wc$ymag,itype="two-ways")
        if(wt_opt==TRUE){
          a <- opt_wt(xi,yi,wmag,max_wt_exp=wt_opt_max_exp,steps=wt_opt_steps,order=1)
          opt_pVal <- as.double(a$pVal)
          wt_exp <- a$wt_exp
        }
        vxi <- as.vector(xi)
        vyi <- as.vector(yi)
        fit=lm(na.omit(vyi)~poly(na.omit(vxi),1),weight=wmag^wt_exp)
        pVal <- as.double(anova(fit)$'Pr(>F)'[1])
        rsquare <- summary(fit)$r.squared
        rsquare_lci <- as.double(CIr(r=rsquare, n = length(na.omit(as.zoo(xi))), level = .90))[1]
        rg <- range(na.omit(index(na.omit(wc$x))))[2]-range(na.omit(index(na.omit(wc$x))))[1]
        n = length(na.omit(as.zoo(xi)))
        ddx_new <- ts.fusion.correlation_A3(ddx,ddy,itype=itype,wt_type=wt_type,wt_exp=wt_exp,wt_opt=FALSE,wt_opt_steps=wt_opt_steps,wt_opt_max_exp=wt_opt_max_exp,maxdt=maxdt,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)              
      } else {
        rsquare <- NA
        pVal <- NA
        opt_pVal <- NA
        rsquare_lci <- NA
        n <-NA
        wt_exp <- NA
        rg <- NA
        #QF == 2 if ts returns less than 2 correlation pairs
        QF <- 2
        ddx_new <- ddx
        ddx_new[] <- NA
      }
      
    } else {
      rsquare <- NA
      pVal <- NA
      opt_pVal <- NA
      rsquare_lci <- NA
      n <-NA
      wt_exp <- NA
      rg <- NA
      #QF == 2 if ts returns less than 2 correlation pairs
      QF <- 2
      ddx_new <- ddx
      ddx_new[] <- NA
    }
    return(c(QF,rg,n,rsquare,rsquare_lci,wt_exp,pVal,ddx_new))
    #return(sum(na.omit(y)))
  }
  
  
  funx <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    length_overlap <- length(bfastts(as.vector(rx[1]),dx,type=c("irregular")))+7
    res <- array(NA,c(length(i),length_overlap))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      # res[i, ] <- t(apply(y[i, ], l, calc_c))
      res[i, ] <- t(apply(y[i, ], 1,get_x))        
    }
    return(res)
  }
  
  b <- brick(addLayer(rx,ry))
  if(cores==1){
    x <- calc(b, fun=funx)
  }else {
    x <- mc.calc(b, fun=funx, mc.cores=cores)
  }
  #names(x) <- index(bfastts(as.vector(rx[1]),dx,type=c("irregular")))
  names(x)[1:7] <- c("quality flag","overlap years","n","r-square","LCIr-square90","wt_exp","p-value")
  
  return(x)
}

get.brick_ts.fusion.correlation_flags_A3_dates <- function(rx,ry,dx,dy,itype=c("one-way","two-ways"),wt_type = c("none", "magnitude"),wt_exp=0,wt_opt=FALSE,wt_opt_steps=1,wt_opt_max_exp=5,maxdt = NULL,interpolation_x=c('linear','season'),interpolation_y=c('linear','season'),order_x=1,order_y=1,maxbreaks_x = 2, maxbreaks_y = 2,cores=1){
  #getfusion brick with flags, wt_exp gets calculated and comitted to ts.fusion.correlation_A3
  #therfore wt_opt is set to always FALSE!
  print(paste("itype: ",itype))
  print(paste("wt_type: ",wt_type))
  print(paste("wt_exp: ",wt_exp))
  print(paste("wt_opt: ",wt_opt))
  print(paste("wt_opt_steps: ",wt_opt_steps))
  print(paste("wt_opt_max_exp: ",wt_opt_max_exp))
  #get x ts for brick cell
  
  ts_to_Date_leapyears <- function(x)
    #bfastts does not consider leapyear!!
    #unique(1900 + as.POSIXlt(dates)$year + (yday365(dates) - 1)/365)
    #length(1900 + as.POSIXlt(dates)$year + (yday365(dates) - 1)/365)
    #duplicated(1900 + as.POSIXlt(dates)$year + (yday365(dates) - 1)/365)
    #find duplicates in ts with refere to 29th February (leapyear)
    #instead of 29th feb bfastts writes 1.march
  {
    #if lapyera add 1
    yu<-function(yr){
      u <- 0
      if(yr %% 4 == 0) {u <- 1}
      return(u)
    }
    conv.frac.yr <- function(yr) as.Date((yr-1970)*(365 + yu(yr)) + round((yr-1970)/4) , origin="1970-01-01" ) 
    return(as.Date(sapply(x, conv.frac.yr), origin="1970-01-01") )
  }
  
  
  get_x <- function(v){
    #split up combined input ts
    #split up combined input ts
    x <- v[1:length(dx)]
    y <- v[(length(dx)+1):(length(dx)+length(v))]
    #create ts
    ddx <- bfastts(x,dx,type=c("irregular"))
    ddy <- bfastts(y,dy,type=c("irregular"))
    #get correlation pairs 
    #cb <- get.ts.correlation(ddx,ddy-ddz,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
    #Check if (i) TS overlap (ii) TS at least 2 valid data points
    
    QF <- 0    
    if(min(ts.overlap_nr(ddx,ddy)[]) > 1){
      wc <- ts.correlation_A3(ddx,ddy,maxdt=maxdt,itype=itype,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      xi <- wc$xi
      yi <- wc$yi
      #ddx_new <- ts.fusion.correlation_A3(ddx,ddy,itype=itype,wt_type=wt_type,wt_exp=wt_exp,maxdt=maxdt,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)
      if(length(na.omit(as.zoo(xi))) >= 2){
        wmag <- get.magweight(wc$x,wc$y,wc$xmag,wc$ymag,itype="two-ways")
        if(wt_opt==TRUE){
          a <- opt_wt(xi,yi,wmag,max_wt_exp=wt_opt_max_exp,steps=wt_opt_steps,order=1)
          #this is a work around and should be fixed in opt_wt
          if(length(a$wt_exp)!=0){
            opt_pVal <- as.double(a$pVal)
            wt_exp <- a$wt_exp
          } else {
            wt_exp <- 0
            opt_pVal <- 1
          }
          
          
        }
        print(wt_exp)
        vxi <- as.vector(xi)
        vyi <- as.vector(yi)
        fit=lm(na.omit(vyi)~poly(na.omit(vxi),1),weight=wmag^wt_exp)
        pVal <- as.double(anova(fit)$'Pr(>F)'[1])
        rsquare <- summary(fit)$r.squared
        rsquare_lci <- as.double(CIr(r=rsquare, n = length(na.omit(as.zoo(xi))), level = .90))[1]
        rg <- range(na.omit(index(na.omit(wc$x))))[2]-range(na.omit(index(na.omit(wc$x))))[1]
        n = length(na.omit(as.zoo(xi)))
        ddx_new <- ts.fusion.correlation_A3(ddx,ddy,itype=itype,wt_type=wt_type,wt_exp=wt_exp,wt_opt=FALSE,wt_opt_steps=wt_opt_steps,wt_opt_max_exp=wt_opt_max_exp,maxdt=maxdt,interpolation_x=interpolation_x,interpolation_y=interpolation_y,order_x=order_x,order_y=order_y,maxbreaks_x=maxbreaks_x,maxbreaks_y=maxbreaks_y)        
        
        #work around to decrease ts to actual observation of x and y ts
        #combined ts
        tsi <- as.vector(c(rx[1],ry[1]))
        tsi[] <- 1
        #combined dates
        fus_dates <- c(dx, dy)
        tsi_bfastts <- bfastts(tsi,fus_dates,type=c("irregular"))
        #reduce ts to elements for which x or y obs exisited
        index <- which(tsi_bfastts == 1)
        ddx_new <- ddx_new[index]       
        
      } else {
        rsquare <- NA
        pVal <- NA
        opt_pVal <- NA
        rsquare_lci <- NA
        n <-NA
        wt_exp <- NA
        rg <- NA
        #QF == 2 if ts returns less than 2 correlation pairs
        QF <- 2
        #ddx_new <- ddx
        #ddx_new[] <- NA
        ddx_new <- rep(NA, length(as.vector(c(rx[1],ry[1]))))
      }
      
    } else {
      rsquare <- NA
      pVal <- NA
      opt_pVal <- NA
      rsquare_lci <- NA
      n <-NA
      wt_exp <- NA
      rg <- NA
      #QF == 2 if ts returns less than 2 correlation pairs
      QF <- 2
      #ddx_new <- ddx
      #ddx_new[] <- NA
      ddx_new <- rep(NA, length(as.vector(c(rx[1],ry[1]))))
    }
    print(c(QF,rg,n,rsquare,rsquare_lci,wt_exp,pVal,length(ddx_new)))
    return(c(QF,rg,n,rsquare,rsquare_lci,wt_exp,pVal,ddx_new))
    #return(sum(na.omit(y)))
  }
  
  
  funx <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    length_overlap <- length(as.vector(c(rx[1],ry[1])))+7
    res <- array(NA,c(length(i),length_overlap))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      # res[i, ] <- t(apply(y[i, ], l, calc_c))
      res[i, ] <- t(apply(y[i, ], 1,get_x))        
    }
    return(res)
  }
  
  b <- brick(addLayer(rx,ry))
  if(cores==1){
    x <- calc(b, fun=funx)
  }else {
    x <- mc.calc(b, fun=funx, mc.cores=cores)
  }
  #names(x) <- index(bfastts(as.vector(rx[1]),dx,type=c("irregular")))
  names(x)[1:7] <- c("quality flag","overlap years","n","r-square","LCIr-square90","wt_exp","p-value")
  
  #work around to assign Landsat names
  tsi <- as.vector(c(rx[1],ry[1]))
  tsi[] <- 1
  fus_dates <- c(dx, dy)
  tsi_bfastts <- bfastts(tsi,fus_dates,type=c("irregular"))
  index <- which(tsi_bfastts == 1)
  #dx_bfastts <- index(bfastts(as.vector(rx[1]),dx,type=c("irregular")))
  #dx_bfastts <- dx_bfastts[index]
  tsi_bfastts_dates <- index(tsi_bfastts)
  tsi_bfastts_dates <- tsi_bfastts_dates[index]
  
  names(x)[8:nlayers(x)] <- paste("LE7075072",format(ts_to_Date_leapyears(tsi_bfastts_dates),"%Y"),format(ts_to_Date_leapyears(tsi_bfastts_dates),"%j"),sep="")
  return(x)
}

