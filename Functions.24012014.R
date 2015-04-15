require(bfast) 

#Pixel based functions (copied from fusion work correlation)

##############################################
#Plot functions
##############################################

#Function plot.bfastmonitor_new from BFAST package
plot.bfastmonitor_new <- function(x, na.rm = TRUE, main = TRUE, ylab = "Data", ...){
  if(isTRUE(main)) main <- if(is.na(x$breakpoint)) {
    "No break detected"
  } else {
    sprintf("Break detected at: %i(%i)", floor(x$breakpoint),
            round((x$breakpoint - floor(x$breakpoint)) * frequency(x$data)) + 1)
  }
  
  y <- if(is.null(dim(x$data))) x$data else x$data[,1L]
  if(na.rm) y <- na.omit(as.zoo(y))
  plot(y, type = "n", main = main, ylab = ylab, ...)
  lines(window(y, end = x$history[2]), col = "black")
  lines(window(y, start = x$history[1], end = x$history[2]),
        col = "darkgreen", type = "p", pch = 19, cex = 0.5)
  lines(window(y, start = x$monitor[1]), col = "red")
  points(window(y, start = x$monitor[1]), col = "red", pch=19, cex=0.5) # new
  test_pred <- predict(x$model, newdata = x$tspp)
  test_pred <- zoo(test_pred, x$tspp$time, frequency = frequency(y))
  lines(test_pred, col = "blue", lwd = 1.5)
  
  abline(v = x$monitor[1], lty = 2, col = "black", lwd = 1)
  abline(v = x$breakpoint, lty = 2, col = "red", lwd = 2)
  
  legend("bottomleft", bty = "n",
         c("Historical data", "New data", "Stable history", "Fit based on stable history", "Start of the Monitoring period", "Time of detected break"),
         lty = c(1, 1, NA, 1, 2, 2),
         col = c("black", "red", "darkgreen", "blue", "black", "red"),
         pch = c(NA, NA, 19, NA, NA, NA)
  )
  invisible(x)
}
#Function - Plot PALSAR single pixel TS (HH, HV and HHHVratio)
plot.PALSAR_single <- function(hh,hv,title="",mf = c(1,1)) {
  #hh and hv are expected to be irregular ts generated with bfastts
  
  #omit NAs from irregular TS
  hh <- na.omit(as.zoo(hh))
  hv <- na.omit(as.zoo(hv))
  
  #calculate HHHV ratio
  hhhv = hh - hv
  opar <- par(mfrow=mf)
  
  #HH and HV plot
  p <- plot(hh, type='b', col="red",ylim=c(-30,-5), ylab="dB", xlab="Acquisition date", main = title, cex.main = 1)
  p <- p + lines(hv, type='b', col="green")
  legend("bottomleft", c(paste("HH m:",round(mean(hh),2)," sd:",round(sd(hh),2)), 
                         paste("HV m:",round(mean(hv),2)," sd:",round(sd(hv),2))), 
         cex=0.8, col=c("red","green"), pch=21:21, lty=1:1);
  
  #subplot HHHV ratio
  empty <- zoo(rep(0,19),time(hh)) # empty dataset with HH time stamp
  plot(empty, ylim=c(-5,15),xlab="Acquisition date",ylab="dB")
  lines(hhhv, type='b') 
  legend("bottomleft", cex=0.8, paste("HHHVratio m:",round(mean(hhhv),2)," sd:",round(sd(hhhv),2)), pch=21:21, lty=1:1);
  
  #set par back to 1/1
  par(mfrow=c(1,1))
}
#Function to plot HHHVratio on top of NDVI
plot.NDVI_HH_HV <- function(u,bfm_u1,title="") {
  #omit NAs from irregular TS
  ndvi <- na.omit(as.zoo(u[,1]))
  hh <- na.omit(as.zoo(u[,2]))
  hv <- na.omit(as.zoo(u[,3]))
  ndvi_obs <- na.omit(as.zoo(u[,5]))
  #cloud cover
  cc <- 1-length(na.remove(ndvi))/length(na.remove(ndvi_obs))
  #set par
  par(mar=c(5, 4, 4, 6) )
  
  # !set y data range from data range!
  ylim_1 = c(range(ndvi))
  ylim_2 = c(range(hh))
  ylim_3 = c(range(hv))
  # !set y data range manually
  #ylim_1 = c(0,1)
  #ylim_2 = c(-18,-5)
  #ylim_3 = c(-26,-10)
  
  #slope per step *frequency = slope per year
  slp_yr <- bfm_u1$model$coefficients[2]*frequency(time(bfm_u1$data))
  
  #plot NDVI
  p <- plot(ndvi,xlim=range(time(u)),axes=F, ylim=ylim_1, xlab="", ylab="",type="l",col="black",pch=".")
  p <- p + points(ndvi,xlim=range(time(u)),ylim=ylim_1, pch = 19,cex=0.5)
  p <- axis(2, ylim=ylim_1,col="black",lwd=1)
  mtext(2,text="NDVI",line=2)
  axis(1,pretty(range(time(u)),10))
  mtext(3,text=paste("ID: ", title," ;Cloud cover: ",round(cc, digits=2),"; Slope/year :",round(slp_yr, digits=4),sep=""),line=1)
  #plot trend according to plot.bfastmonitor
  y <- if(is.null(dim(bfm_u1$data))) bfm_u1$data else bfm_u1$data[,1L]
  test_pred <- predict(bfm_u1$model, newdata = bfm_u1$tspp)
  test_pred <- zoo(test_pred, bfm_u1$tspp$time, frequency = frequency(y))
  lines(test_pred, col = "black", lwd = 1.5)
  #plot HH
  par(new=T)
  p <- plot(hh,xlim=range(time(u)),ylim=ylim_2, axes=F, xlab="", ylab="",type="l",col="blue",pch=".")
  p <- p + points(hh,xlim=range(time(u)),ylim=ylim_2, pch = 19,cex=0.5,col="blue")
  p <- axis(4, ylim=ylim_2,col="blue",lwd=1)
  mtext(4,text="HH",line=2,col="blue")
  #plot HV
  par(new=T)
  p <- plot(hv,xlim=range(time(u)),ylim=ylim_3, axes=F, xlab="", ylab="", type="l",col="red",pch=".")
  p <- p + points(hv,xlim=range(time(u)),ylim=ylim_3, pch = 19,cex=0.5,col="red")
  p <- axis(4, ylim=ylim_3,col="red",lwd=1,line=3)
  mtext(4,text="HV",line=5,col="red")
}
#Function to plot HHHVratio on top of NDVI
plot.NDVI_HHHVratio <- function(u,bfm_u1,title=""){
  #omit NAs from irregular TS
  ndvi <- na.omit(as.zoo(u[,1]))
  hhhv <- na.omit(as.zoo(u[,4]))
  ndvi_obs <- na.omit(as.zoo(u[,5]))
  #cloud cover
  cc <- 1-length(na.remove(ndvi))/length(na.remove(ndvi_obs))
  #set par
  par(mar=c(5, 4, 4, 6) )
  # !set y data range from data range!
  par(mar=c(5, 4, 4, 6))
  ylim_1 = c(range(ndvi))
  ylim_4 = c(range(hhhv*-1))
  #ylim_1 = c(0,1)
  #ylim_2 = c(-10,-4)
  
  #slope per step *frequency = slope per year
  slp_yr <- bfm_u1$model$coefficients[2]*frequency(time(bfm_u1$data))
  
  
  #plot NDVI
  p <- plot(ndvi,xlim=range(time(u)),axes=F, ylim=ylim_1, xlab="", ylab="",type="l",col="black",pch=".")
  p <- p + points(ndvi,xlim=range(time(u)),ylim=ylim_1, pch = 19,cex=0.5)
  p <- axis(2, ylim=ylim_1,col="black",lwd=1)
  mtext(2,text="NDVI",line=2)
  axis(1,pretty(range(time(u)),10))
  mtext(3,text=paste("ID: ", title," ;Cloud cover: ",round(cc, digits=2),"; Slope/year :",round(slp_yr, digits=4),sep=""),line=1)
  #plot trend according to plot.bfastmonitor
  y <- if(is.null(dim(bfm_u1$data))) bfm_u1$data else bfm_u1$data[,1L]
  test_pred <- predict(bfm_u1$model, newdata = bfm_u1$tspp)
  test_pred <- zoo(test_pred, bfm_u1$tspp$time, frequency = frequency(y))
  lines(test_pred, col = "black", lwd = 1.5)
  #plot HHHVratio
  par(new=T)
  p <- plot(hhhv*-1,xlim=range(time(u)),ylim=ylim_4, axes=F, xlab="", ylab="",type="l",col="blue",pch=".")
  p <- p + points(hhhv*-1,xlim=range(time(u)),ylim=ylim_4, pch = 19,cex=0.5,col="blue")
  p <- axis(4, ylim=ylim_4,col="blue",lwd=1)
  mtext(4,text="HHHVratio*-1",line=2,col="blue")
}
#Function to plot two time-series
plot.2ts <- function(x,y,ylim_x=NULL,ylim_y=NULL,name_x="",name_y="",col_x="black",col_y="blue", points=TRUE, title="") {
  #omit NAs from irregular TS
  zx <- na.omit(as.zoo(x))
  zy <- na.omit(as.zoo(y))
  #set par
  par(mar=c(5, 4, 4, 6) )
  
  # !set y data range from data range!
  if (is.null(ylim_x)){
    ylim_x = c(range(zx))
    #ylim_y = c(range(zy))
  } else {
    ylim_x <- ylim_x
  }
  if (is.null(ylim_y)){
    ylim_y = c(range(zy))
    #ylim_y = c(range(zy))
  } else {
    ylim_y <- ylim_y
  }
  
  #plot x
  p <- plot(zx,xlim=range(time(x)),axes=F, ylim=ylim_x, xlab="", ylab="",type="l",col=col_x,pch=".")
  if(points==TRUE) {p <- p + points(zx,xlim=range(time(x)),ylim=ylim_x, pch = 19,cex=0.5)}
  p <- axis(2, ylim=ylim_x,col=col_x,lwd=1)
  mtext(2,text=name_x,line=2)
  axis(1,pretty(range(time(x)),10))
  mtext(3,text=paste(title,sep=""),line=1)
  #plot y
  #cat ("Press [enter] to plot 2nd time-series")
  #line <- readline()
  par(new=T)
  p <- plot(zy,xlim=range(time(x)),ylim=ylim_y, axes=F, xlab="", ylab="",type="l",col=col_y,pch=".")
  if(points==TRUE) {p <- p + points(zy,xlim=range(time(x)),ylim=ylim_y, pch = 19,cex=0.5,col="blue")}
  p <- axis(4, ylim=ylim_y,col=col_y,lwd=1)
  mtext(4,text=name_y,line=2,col=col_y)
}
#Function to plot two time-series
plot.1ts <- function(x,ylim_x=NULL,name_x="",points=TRUE, title="",col_x="black") {
  #omit NAs from irregular TS
  zx <- na.omit(as.zoo(x))
  #set par
  par(mar=c(5, 4, 4, 6) )
  
  # !set y data range from data range!
  if (is.null(ylim_x)){
    ylim_x = c(range(zx))
    #ylim_y = c(range(zy))
  } else {
    ylim_x <- ylim_x
  }
  
  #plot x
  p <- plot(zx,xlim=range(time(x)),axes=F, ylim=as.double(ylim_x), xlab="", ylab="",type="l",col=col_x,pch=".")
  if(points==TRUE) {p <- p + points(zx,xlim=range(time(x)),ylim=as.double(ylim_x), pch = 19,cex=0.5,col=col_x)}
  p <- axis(2, ylim=ylim_x,col=col_x,lwd=1)
  mtext(2,text=name_x,line=2)
  axis(1,pretty(range(time(x)),10))
  mtext(3,text=paste(title,sep=""),line=1)
  #plot y
  #cat ("Press [enter] to plot 2nd time-series")
  #line <- readline()
  
}

##############################################
#TS Brick functions
##############################################
#get brick with
get.brick.tsMD_random <- function(x,percMD,cores=1){
  #get x ts for brick cell
  get_x <- function(v){
    vcc <- get.tsMD_random(v,percMD)
    #reduce time-series density to xx no-data
    return(vcc)
    #return(sum(na.omit(y)))
  }
  
  funx <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    length_overlap <- nlayers(x)
    res <- array(NA,c(length(i),length_overlap))
    # do not apply if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      # res[i, ] <- t(apply(y[i, ], l, calc_c))
      res[i, ] <- t(apply(y[i, ], 1,get_x))        
    }
    return(res)
  }
  
  
  if(cores==1){
    r <- calc(x, fun=funx)
  }else {
    r <- mc.calc(x, fun=funx,mc.cores=cores)
  }
  return(r)
}
#Inverse spatial buffer
inverse_buffer <- function(y,range=c(1990,2015), width=30){
  #new raster
  ny <- y
  ny[!is.na(ny)] <- NA
  for(i in range[1]:range[2]) {
    #work raster
    wy <- y
    #set all raster values != i to NA
    wy[y!=i] <- NA
    #check wether i value exisit in raster
    v <- as.vector(wy)
    if(length(v[!is.na(v)])>0){
      #invert NA and raster value
      wy[is.na(wy)] <- 1
      wy[wy!=1]<-NA
      #buffer actual NA vaules 
      wy <- buffer(wy,width=width)
      #assign
      ny[is.na(wy)] <- i
    }
  }
  return(ny)
}
#function that returns a raster layer with the correlation for each cell
#Calculates the percentage of NA per pixel for a raster brick
calc.brick.percNA <- function(x,cores=1){
  fun <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    return(percNA)
  }
  
  if(cores==1){
    r <- calc(x, fun=fun)
  }else {
    r <- mc.calc(x, fun=fun,mc.cores=cores)
  }
  return(r)
}
calc.brick.perc99 <- function(x,cores=1){
  fun <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(x==-99)/length(x)))
    return(percNA)
  }
  
  if(cores==1){
    r <- calc(x, fun=fun)
  }else {
    r <- mc.calc(x, fun=fun,mc.cores=cores)
  }
  return(r)
}
calc.brick.percInf <- function(x,cores=1){
  fun <- function(y){
    percNA <- apply(y, 1, FUN = function(x) (sum(x==Inf)/length(x)))
    return(percNA)
  }
  
  if(cores==1){
    r <- calc(x, fun=fun)
  }else {
    r <- mc.calc(x, fun=fun,mc.cores=cores)
  }
  return(r)
}



##############################################
#TS functions
##############################################
#Check if valid observation (not NA) in 
#two TS are overlapping, return TRUE if overlapping
ts.overlap <- function(x,y){
  zox <- na.omit(as.zoo(x))
  zoy <- na.omit(as.zoo(y))
  if (end(zox) > start(zoy) && end(zoy) > start(zox)) {

    ov <- TRUE
  } else {
    ov <- FALSE
  }
  return(ov)
}
ts.overlap_nr <- function(x,y){
  zox <- na.omit(as.zoo(x))
  zoy <- na.omit(as.zoo(y))
  if (end(zox) > start(zoy) && end(zoy) > start(zox)) {
    
    ov <- TRUE
    #concert days to zoo
    zx <- as.zoo(x)
    zy <- as.zoo(y)
    #start and end of overlappping period
    st <- max(start(na.omit(zx)),start(na.omit(zy)))
    en <- min(end(na.omit(zx)),end(na.omit(zy)))
    #include one valid observation exceeding the overlapping ts by maxdt days (max and min) 
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
    return(c(length(na.omit(zxw)),length(na.omit(zyw))))
  } else {
    return(c(0,0))
  }
  
}
#Calculates the percentage of NA per pixel for a vector
calc.percNA <- function(x,cores=1){
  return((sum(is.na(x))/length(x)))
}
#Function for simulate cloud cover and reduce observations
cloud_cover <- function(ts,ts_obs,cc_i){
  #number of valid observations
  l_ts <-  length(na.remove(ts))
  #number of all observations
  l_ts_obs <- length(na.remove(ts_obs))
  #calculate how x-te valid observation must be removed to meet wanted cloud cover
  ag <- l_ts/(l_ts - l_ts_obs*(1-cc_i))
  print(ag)
  ts_cc <- ts
  #set current step in valid observations
  a <- 1
  for(i in 1:length(ts)) {
    #check if valid observatoin
    if(!is.na(ts_cc[i])){
      #if step > x-te obs to be removed
      if (a>=ag){
        #print(a)
        #remove observation
        ts_cc[i] <- NA
        #increase with 1 step and decrease x-te obs that was just removed
        a <- 1 + (a-ag)
        #print(a)
      } else {
        a <- a + 1
      } 
    } 
  }
  return(ts_cc)
}
#Function for increasing percentage of missing data (MD) in ts by randomly excluding valid observation
#If percMD < exising percentageMD than the original ts gets returned
get.tsMD_random <- function(ts,percMD){
  #number of ts observations
  ts_obs <- ts  
  ts_obs[] <- 1
  #number of valid observations
  l_ts <-  length(na.omit(ts))
  #number of all observations
  l_ts_obs <- length(na.omit(ts_obs))
  #number of observations to be removed
  ag <- l_ts - l_ts_obs*(1-percMD)
  #tsMD = ts with missing data
  tsMD <- ts
  #if ag > 0 (existing MD < cc_i)
  if(ag > 0){
    #Select ag (number of valid observations to be removed) random numbers between 1 and number of valid observations, without replacement
    x <- sample(1:l_ts, ag, replace=F)
    for(i in 1:length(x)) {
      #which(index(ts)==index(na.omit(as.zoo(ts)))[x[i]])
      tsMD[which(index(as.zoo(ts))==index(na.omit(as.zoo(ts)))[x[i]])] <- NA
    } 
  }   
  return(tsMD)
}
#calc.Rsquare
calc.Rsquare.plot <- function(x,y,order=1,name_x="",name_y=""){
  p <- plot(x,y,xlab=name_x,ylab=name_y)
  vx <- as.vector(x)
  vy <- as.vector(y)
  fit=lm(na.omit(vy)~poly(na.omit(vx),order))
  r <- range(na.omit(vx))
  xx <- seq(r[1],r[2], length.out=250)
  predict(fit, data.frame(vx=xx))
  lines(xx, predict(fit, data.frame(vx=xx)), col='blue')
  title(paste("R²: ", round(summary(fit)$r.squared,digits=4),"; Adj. R²: ",round(summary(fit)$adj.r.squared,digits=4),"; p-value: ",formatC(anova(fit)$'Pr(>F)'[1],digits=12,format="f")," [",length(na.omit(x)),"]",sep=""),cex.main=1)
  return(list(r.squared = summary(fit)$r.squared,adj.r.squared = summary(fit)$adj.r.squared))
}
calc.Rsquare <- function(x,y,order=1){
  vx <- as.vector(x)
  vy <- as.vector(y)
  fit=lm(na.omit(vy)~poly(na.omit(vx),order))
  return(list(r.squared = summary(fit)$r.squared,adj.r.squared = summary(fit)$adj.r.squared))
}
#calc.Rsquare with weight
calc.Rsquare_weight.plot <- function(x,y,weight,order=1,name_x="",name_y="",ylim=NULL,xlim=NULL){
  p <- plot(x,y,xlab=name_x,ylab=name_y,xlim=xlim,ylim=ylim)
  vx <- as.vector(x)
  vy <- as.vector(y)
  fit=lm(na.omit(vy)~poly(na.omit(vx),order),weight=weight)
  r <- range(na.omit(vx))
  xx <- seq(r[1],r[2], length.out=250)
  predict(fit, data.frame(vx=xx))
  lines(xx, predict(fit, data.frame(vx=xx)), col='blue')
  title(paste("r²: ", round(summary(fit)$r.squared,digits=4),";   p-value: ",formatC(anova(fit)$'Pr(>F)'[1],digits=10,format="f"),"   [n=",length(na.omit(x)),"]",sep=""),cex.main=1)
  return(list(r.squared = summary(fit)$r.squared,adj.r.squared = summary(fit)$adj.r.squared))
}
calc.Rsquare_weight <- function(x,y,weight,order=1){
  vx <- as.vector(x)
  vy <- as.vector(y)
  fit=lm(na.omit(vy)~poly(na.omit(vx),order),weight=weight)
  return(list(r.squared = summary(fit)$r.squared,adj.r.squared = summary(fit)$adj.r.squared))
}


calc.Slope <- function(x,y,order=1){
  vx <- as.vector(x)
  vy <- as.vector(y)
  fit=lm(na.omit(vy)~poly(na.omit(vx),order))
  return(summary(fit)$coef[2])
}
calc.Slope_weight <- function(x,y,weight,order=1){
  vx <- as.vector(x)
  vy <- as.vector(y)
  fit=lm(na.omit(vy)~poly(na.omit(vx),order),weight=weight)
  return(summary(fit)$coef[2])
}

##############################################
#Backcast functions
##############################################
#Backcast SAR TS for dates for which NDVI values exisit
ts.backcast <- function(u,start_overlap,end_overlap,bfm){
    uo <- window(u, start = start_overlap, end = end_overlap)
    #for-loop that checks if [u,1] has a value and if <2007 (<overlapping period do not to overwrite values from overlapping perio)
    #and if so it calcualtes values for HV -> mean value of 2007 +- random value around sd of 2007
    bc <- rep(2)
    bc[1] <- bfm$history[1]
    bc[2] <- max(time(window(u[,1], end = start_overlap)))
    for(i in 1:length(u[,1])) {
      if(!is.na(u[i,1])&time(u[,1])[i]>=bc[1] & time(u[,1])[i]<=bc[2]){
        u[i,2] <- runif(1,mean(uo[,2],na.rm=TRUE)-sd(uo[,2],na.rm=TRUE),mean(uo[,2],na.rm=TRUE)+sd(uo[,2],na.rm=TRUE))
        u[i,3] <- runif(1,mean(uo[,3],na.rm=TRUE)-sd(uo[,3],na.rm=TRUE),mean(uo[,3],na.rm=TRUE)+sd(uo[,3],na.rm=TRUE))
        u[i,4] <- runif(1,mean(uo[,4],na.rm=TRUE)-sd(uo[,4],na.rm=TRUE),mean(uo[,4],na.rm=TRUE)+sd(uo[,4],na.rm=TRUE))
      }
    }
    return(u)
  }
  
#Interpolate historical part of TS if historical part is too short
#Normal interpolation function without considering the sd
ts.fill <- function(ts,num_fill,start_over){
    #end of ts part for which data is interpolated
    ts_ap_end <- max(time(window(ts, end = start_over)))
    #interpolate data between min and max timestep time series part to be reconstructed
    #min and max time steps constitute the min and max border
    #num_fill = final number of time-steps after interpolation
    ap <- approx(ts[which(time(ts)[]<=ts_ap_end )],n=num_fill)
    #approx returns a list. unlist returns a vector based on a list
    ap_x_v <- unlist(ap[1]) #index in ts
    ap_y_v <- unlist(ap[2]) #data
    ts[ap_x_v] <- ap_y_v
    return(ts)
  }
#Interpolates not just in betwen the ts, but in the range of mean-sd and mean+sd
#Crucial to maintain the sd characteristics. Otherwise, BFAST has problems to find a stable history
ts.fill_sd <- function(ts,num_fill,start_over,end_over){
    #end of ts part for which data is interpolated
    ts_ap_end <- max(time(window(ts, end = start_over)))
    #interpolate data between min and max timestep time series part to be reconstructed
    #min and max time steps constitute the min and max border
    #num_fill = final number of time-steps after interpolation
    ap <- approx(ts[which(time(ts)[]<=ts_ap_end)],n=num_fill)
    #approx returns a list. unlist returns a vector based on a list
    ap_x_v <- unlist(ap[1])
    ap_y_v <- unlist(ap[2])
    #Calculate sd and mean of overlapping period in TS
    ts_over <- window(ts, start = start_overlap, end = end_overlap)
    ts_over_sd <- sd(ts_over,na.rm=TRUE)
    ts_over_mean <- mean(ts_over,na.rm=TRUE) 
    #Interpolates not just in betwen the ts, but in the range of mean-sd and mean+sd
    ts[ap_x_v] <- unlist(lapply(ap_y_v,function(x) runif(1,ts_over_mean-ts_over_sd,ts_over_mean+ts_over_sd)))
    #ts[ap_x_v] <- ap <- runif(1,,na.rm=TRUE)-sd(u07,na.rm=TRUE),mean(u07,na.rm=TRUE)+sd(u07,na.rm=TRUE))
    return(ts)
  }
  

##################################################
# bfastmonitor reich006 function
#################################################

bfastmonitor_reich006 <- function(data, start,
                                  formula = response ~ trend + harmon,
                                  order = 3, lag = NULL, slag = NULL,
                                  history = c("ROC", "BP", "all"),
                                  type = "OLS-MOSUM", h = 0.25, end = 10, level = 0.05,
                                  hpc = "none", verbose = FALSE, plot = FALSE, hdyn = FALSE,history_obs_min=2)
  ##reich006: Adapted bfastmonitor to run it on spare ts with history period >= 2
  ##reich006: Change 1: line27: Check if at least 1 obs in history period, otherwise stop function and return NA
  ##reich006: Change 2: line60: Dynamic h-value (for h=0.25 at least 8 observations in historic period have to exist, h=0.5 for at least 4)
  ##reich006: !!! this could also be h=0.25 for >=4 and h=0.5 >=2 if floor(h * NROW(test_tspp)) < 1 instead of <=1!
  ##reich006: !!! Discuss with JV
  ##reich006: Change 3: line75: check for minimum number of observations in history; >= 2 seem to be threshold that method provides change
{
  ## PREPROCESSING
  ## two levels needed: 1. monitoring, 2. in ROC (if selected)
  level <- rep(level, length.out = 2)
  
  if(!is.ts(data)) data <- as.ts(data)
  
  ## frequency of data
  freq <- frequency(data)
  ## start on natural scale (if necessary)
  time2num <- function(x) if(length(x) > 1L) x[1L] + (x[2L] - 1)/freq else x
  start <- time2num(start)
  
  ## full data
  data_tspp <- bfastpp(data, order = order, lag = lag, slag = slag)
  
  ## reich006: check if at least 1 obs in history, otherwise stop function and return NA 
  if (start <= min(data_tspp[]$time)) {
    warning("Zero observations in history period")
    rval <- list(
      data = data,
      tspp = NA,
      model = NA,
      mefp = NA,
      history = NA,
      monitor = NA,
      breakpoint = NA,
      magnitude = NA
    )
    class(rval) <- "bfastmonitor"
    return(rval)
  }
  
  ## SELECT STABLE HISTORY  
  ## full history period
  history_tspp <- subset(data_tspp, time < start)
  
  
  
  ## find start of history period
  ## (may be specified via character, function, or time index directly)
  if(is.null(history)) {
    history <- start(history_tspp$response)
  } else if(all(is.character(history))) {
    history <- match.arg(history)
    history <- switch(history,    
                      "all" = start(history_tspp$response),      
                      "ROC" = history_roc(formula, data = history_tspp, level = level[2]),
                      "BP" = history_break(formula, data = history_tspp, hpc = hpc)
    )
  } else if(all(is.function(history))) {
    history <- history(formula, data = history_tspp)
  }
  history <- time2num(history)
  
  ## compute subset
  history_tspp <- subset(history_tspp, time >= history)
  
  ## output information (if desired)
  if(verbose) {
    cat("\nBFAST monitoring\n\n1. History period\n")
    cat(sprintf("Stable period selected: %i(%i)--%i(%i)\n",
                start(history_tspp$response)[1], start(history_tspp$response)[2],
                end(history_tspp$response)[1], end(history_tspp$response)[2]))
    cat(sprintf("Length (in years): %f\n", NROW(history_tspp)/freq))
  }
  
  ## MODEL HISTORY PERIOD
  test_tspp <- history_tspp
  
  ## reich006: DYNAMIC h-value
  if(isTRUE(hdyn)){
    h = 0.25 #for at least 8 obs in history
    if(length(history_tspp$response) < 8) {h = 0.5} #for at least 4 obs in history
    if(length(history_tspp$response) < 4) {h = 1} #for at least 1 obs in history
    print(paste("h-value: ",h))
  }
  test_mefp <- mefp(formula, data = test_tspp, 
                    type = type, period = end, h = h, alpha = level[1])
  
  # reich006: Check if minimum number of historic observations;
  # reich006: at least 2-3 obs that method provides a change
  test_lm <- lm(formula, data = test_tspp)
  if (length(history_tspp$response) >= history_obs_min) {
    if(floor(h * NROW(test_tspp)) <= 1 | NROW(test_tspp) <= length(coef(test_lm))) {
      ok <- FALSE
      warning("too few observations in selected history period")
    } else {
      ok <- TRUE
    }
  } else {
    ok <- FALSE
    print(paste("too few observations in selected history period: ",length(history_tspp$response)))
    warning("too few observations in selected history period")
  }

    if(verbose) {
      cat("Model fit:\n")
      print(coef(test_lm))
    }

  ## MONITOR CHANGES IN THE MONITORING PERIOD
  test_tspp <- subset(data_tspp, time >= history)
  if(ok) {
    test_mon <- monitor(test_mefp, data = test_tspp, verbose = FALSE)
    tbp <- if(is.na(test_mon$breakpoint)) NA else test_tspp$time[test_mon$breakpoint]
    if(verbose) {
      cat("\n\n2. Monitoring period\n")
      cat(sprintf("Monitoring starts at: %i(%i)\n", floor(start), round((start - floor(start)) * freq) + 1))
      if(is.na(tbp)) {      
        cat("Break detected at: -- (no break)\n\n")
      } else {
        cat(sprintf("Break detected at: %i(%i)\n\n", floor(tbp), round((tbp - floor(tbp)) * freq) + 1))
      }
    }
  } else {
    test_mon <- NA
    test_lm <- NA
    tbp <- NA
  }
  
  ## the magnitude of change
  if(ok) {
    test_tspp$prediction <- predict(test_lm, newdata = test_tspp)
    new_data <- subset(test_tspp, time>=start) ## only data from the monitoring period
    magnitude <- median(new_data$response - new_data$prediction,na.rm=TRUE)
  } else {
    test_tspp$prediction <- NA
    magnitude <- NA
  }
  
  ## set up return object
  rval <- list(
    data = data,
    tspp = test_tspp,
    model = test_lm,
    mefp = test_mon,
    history = c(head(history_tspp$time, 1), tail(history_tspp$time, 1)),
    monitor = c(start, tail(test_tspp$time, 1)),
    breakpoint = tbp,
    magnitude = magnitude
  )
  class(rval) <- "bfastmonitor"
  
  ## plot if desired
  if(plot) plot(rval)
  
  ## return object
  return(rval)
}

run_bfm_reich006_cloud_mc <- function(data, history="ROC", monperiod=c(), monend="full", formula=response~trend+harmon, order=1, type = "OLS-MOSUM", 
                                      sceneID=NULL, filename="", overwrite=FALSE,mc.cores=mc.cores, cloud=TRUE,h=0.25,hdyn = FALSE,history_obs_min=2)
  #To run it on SAR TS input needs to be positive (e.g. hv * -1 instead of hv)
{
  
  if(is.null(monperiod)) 
    stop("missing value for monperiod (ie. start of the monitoring period)")
  if(monend[1]!="full" & !is.numeric(monend)) 
    stop("monend should either be \'full\' or numeric of length=2 (last date of monitoring period)")
  #if(monend<monperiod) 
  #stop("monend must be greater than monperiod")
  
  b <- data
  
  # get dates and years from sceneID vector or calculate dates and years from layer names
  if(is.null(sceneID))
    sceneID <- row.names(get.sceneinfo(names(b)))
  dates <- get.sceneinfo(sceneID)$date
  
  # trim time series if monend!="full"
  if(monend[1]!="full"){
    end.date <- as.Date(paste(monend[1], monend[2], sep="-"), format="%Y-%j")
    b <- dropLayer(b, which(dates > end.date))
    
    # redefine dates based on trimmed raster brick
    sceneID <- sceneID[which(dates <= end.date)]
    dates <-get.sceneinfo(sceneID)$date
  }
  
  # declare bfm helper functions
  # TODO: rewrite ybfastmonitor to include bfm.pixel()?
  ybfastmonitor_new <- function(x, dates) {
    vi <- bfastts(x, dates, type = c("irregular"))
    #vi[2] <- vi[1]
    #reich006: vi can also be <= 0 for SAR ts
    #vi[vi <= 0] <- NA
    bfm <- bfastmonitor_reich006(formula = formula, order = order, 
                                 data = vi, start = monperiod, history = history,hdyn=hdyn,history_obs_min=history_obs_min)
    if(cloud==TRUE){
      # if break detected in bfm check if it is valid
      if(!is.na(bfm$breakpoint)){
        vi2 <- vi
        #set break na
        vi2[time(vi2) == bfm$breakpoint] <- NA
        bfm2 <- bfastmonitor_reich006(formula = formula, order = order, 
                                      data = vi2, start = monperiod, history = history,hdyn=hdyn,history_obs_min=history_obs_min)
        #if break detected in bfm2
        if(!is.na(bfm2$breakpoint)){
          #index of detected break bfm1
          bfmi <-  which(time(vi) == bfm$breakpoint)
          #index of detected break bfm2
          bfm2i <- which(time(vi) == bfm2$breakpoint)
          #index of detected break bfm1 in sequence (1,2,3,4 ...) of TS elements with value
          bfm1is <- which(bfmi==which(!is.na(vi)))
          #index of detected break bfm2 in sequence (1,2,3,4 ...) of TS elements with value
          bfm2is <- which(bfm2i==which(!is.na(vi)))
          if(bfm1is+1==bfm2is){
            #print(sprintf("bfm$breakpoint <- bfm$breakpoint"))
            bfm$breakpoint <- bfm$breakpoint
          } else {
            bfm$breakpoint <- bfm2$breakpoint
            #print(sprintf("bfm$breakpoint <- bfm2$breakpoint"))
          }
        } else {
          bfm$breakpoint <- NA
          #print(sprintf("bfm$breakpoint <- NA"))
        }
      }
    }
    #nr of observaitons in history period
    if(length(bfm$history[])>1){lhistory <- length(which((bfm$tspp$time >= bfm$history[1] & bfm$tspp$time <= bfm$history[2]) == TRUE) )
    } else {lhistory <- 0}
    print(lhistory)
    if(!is.na(bfm$breakpoint)){
      break_pred  <- subset(bfm$tspp$prediction,bfm$tspp$time >= bfm$breakpoint)
      break_pred_y <- subset(bfm$tspp$prediction,(bfm$tspp$time >= bfm$breakpoint &  bfm$tspp$time < bfm$breakpoint + 1))
      break_pred
      # 2. get observations of all pixel after breakpoint (and including breakpoint)
      break_response <- subset(bfm$tspp$response,bfm$tspp$time >= bfm$breakpoint)
      break_response_y <- subset(bfm$tspp$response,(bfm$tspp$time >= bfm$breakpoint &  bfm$tspp$time < bfm$breakpoint + 1))
      break_response
      # 3 calc magnitude
      magnitude <- median(break_response - break_pred)
      magnitude_y <- median(break_response_y - break_pred_y)
    } else {
      magnitude <- NA
      magnitude_y <- NA
    }
    
    
    return(cbind(bfm$breakpoint, bfm$magnitude, magnitude, magnitude_y, (bfm$history[2] - 
                                                   bfm$history[1]),lhistory))
  }
  fun <- function(y) {
    percNA <- apply(y, 1, FUN = function(x) (sum(is.na(x))/length(x)))
    i <- ((percNA < 1)) # logical vector corresponding to pixel ts
    res <- matrix(NA, length(i), 6)
    
    # do not apply bfm if there are only NA's in the pixel ts
    if (sum(i) > 0) {
      print("indside")
      res[i, ] <- t(apply(y[i, ], 1, ybfastmonitor_new, dates))
    }
    res
  }
  
  if(filename!="")
    x <- mc.calc(b, fun=fun, filename=filename, overwrite=overwrite, mc.cores=mc.cores)
  else
    x <- mc.calc(b, fun=fun, mc.cores=mc.cores)
  names(x) <- c("bfm.breakpoint","bfm.magnitude","reich006.magnitude","reich006.magnitude_y","bfm.history.years","bfm.history.n")
  return(x)
}

########################################
## Reversely Ordered CUSUM (ROC) test ##
########################################

## A technique to verify whether or not the historical period is stable or not
## reversely order sample and perform
## recursive CUSUM test
history_roc <- function(formula, data, level = 0.05) {
  n <- nrow(data)
  data_rev <- data[n:1,]
  data_rev$response <- ts(data_rev$response)
  y_rcus <- efp(formula, data = data_rev, type = "Rec-CUSUM")
  
  y_start <- if(sctest(y_rcus)$p.value < level) {
    length(y_rcus$process) - min(which(abs(y_rcus$process)[-1] > boundary(y_rcus)[-1])) + 1
  } else {
    1    
  }
  data$time[y_start]
}
