require(raster)
require(zoo)
require(bfast)  ## this loads the bfast package
require(tseries)
#require(PreFAST)
require(caret)
require(psychometric)
require(bfastSpatial)
#require(bayts)

#load all pixel based functions
source("Functions.24012014.R")
source("Functions.correlation.24012014.R")

##################
#load raster data
{
  load("/media//DATA2/reich006/Scripts/R/sar_landsat_ts_fusion_3/data/Data_75_72_2000_2012__750_1400_1000_1650.RData")
#load("data/Data_75_72_2000_2012__750_1400_1000_1650_north.RData")
e_new <- c(563995,572495,8054965,8060465)
#e_new <- c(569995,571495,8054965,8055965)
rndvi <- crop(rndvi,e_new)
rhh <- crop(rhh,e_new)
rhv <- crop(rhv,e_new)
rhvhh <- crop(rhvhh,e_new)
rhh_mt <- crop(rhh_mt,e_new)
rhv_mt <- crop(rhv_mt,e_new)
rhvhh_mt <- crop(rhvhh_mt,e_new)
rndvi60 <- crop(rndvi60,e_new)
rndvi70 <- crop(rndvi70,e_new)
rndvi80 <- crop(rndvi80,e_new)
rndvi90 <- crop(rndvi90,e_new)
rndvi95 <- crop(rndvi95,e_new)
rndvi98 <- crop(rndvi98,e_new)
ls_map <- crop(ls_map,e_new)
ylog <- crop(ylog,e_new)
yplant <- crop(yplant,e_new)
yperiod <- crop(yperiod,e_new)
ync <- crop(ync,e_new)

#restrict to observation period
rndvi <- dropLayer(rndvi,c(1:37))
rndvi <- dropLayer(rndvi,c(81:96))
rndvi95 <- dropLayer(rndvi95,c(1:37))
rndvi95 <- dropLayer(rndvi95,c(81:96))
ndvi_dates <- ndvi_dates[38:133]
ndvi_dates <- ndvi_dates[1:80]
}

#check percentage Missing data in NDVI time series
plot(calc.brick.percNA(rndvi,cores=4))

#Increase missing data pecentage
rndvi90 <- get.brick.tsMD_random(rndvi,0.9,cores=4)
names(rndvi90) <- names(rndvi)
plot(calc.brick.percNA(rndvi90,cores=4))

#Reference data
plot(ylog)


###Getting single pixel
io <- rhv_mt
plot(rndvi,2)
x <- click(io, n=1, id=TRUE, xy=TRUE,cell=TRUE)
#x[3] <- cellFromXY(io,c(565800, 8056630))
#x[1]<-565800
#x[2]<-8056630
hh <- as.vector(rhh[x[1,3]])
hv <- as.vector(rhv[x[1,3]])
hvhh <- as.vector(rhvhh[x[1,3]])
hh_mt <- as.vector(rhh_mt[x[1,3]])
hv_mt <- as.vector(rhv_mt[x[1,3]])
hvhh_mt <- as.vector(rhvhh_mt[x[1,3]])
ndvi <- as.vector(rndvi[x[1,3]])
ndvi_cc <- as.vector(rndvi90[x[1,3]])
dd_hh <- bfastts(hh,hh_dates,type=c("irregular"))
dd_hv <- bfastts(hv,hv_dates,type=c("irregular"))
dd_hvhh <- bfastts(hvhh,hv_dates,type=c("irregular"))
dd_hh_mt <- bfastts(hh_mt,hh_dates,type=c("irregular"))
dd_hv_mt <- bfastts(hv_mt,hv_dates,type=c("irregular"))
dd_hvhh_mt <- bfastts(hvhh_mt,hv_dates,type=c("irregular"))
dd_ndvi <- bfastts(ndvi,ndvi_dates,type=c("irregular"))
dd_ndvi_cc <- bfastts(ndvi_cc,ndvi_dates,type=c("irregular"))

#1) Plot SAR data over NDVI
plot2ts(dd_ndvi,dd_hh_mt,lab_ts1="Landsat NDVI [MD=org]",lab_ts2="PALSAR HH")
plot2ts(dd_ndvi,dd_hv_mt,lab_ts1="Landsat NDVI [MD=org]",lab_ts2="PALSAR HV")
plot2ts(dd_ndvi,dd_hvhh_mt,lab_ts1="Landsat NDVI [MD=org]",lab_ts2="PALSAR HHHV")
plot2ts(dd_ndvi_cc,dd_hvhh_mt,lab_ts1="Landsat NDVI [MD=95]",lab_ts2="PALSAR HHHV")

#################################
#A MulTiFuse - single steps
#################################
xts <- dd_ndvi_cc
yts <- dd_hvhh_mt

#Step 1a: Calculate correlation weight
wc <- calcRegWeight(xts,yts)

#plot original overlapping ts
plot2ts(wc$x,wc$y)
#plot interpolated overlatpping ts
plot2ts(wc$xi,wc$yi)

#Step 1b and c: Regression weight optimization and regression analysis
#opt <- opt_wt.plot(wc$xi,wc$yi,wc$w,max_wt_exp=10,steps=0.1,order=1)
opt <- optimizeRegWeight(wc$xi,wc$yi,wc$w,max_ewf=10,steps=0.1,order=1,plot=TRUE)
#optimized ewf
ewf <- opt$ewf.optimized

weight <- wc$w^0 #no regression weight
weight <- wc$w^1 #simple regression weight
weight <- wc$w^2 #squared regression weight
weight <- wc$w^ewf #optimized regression weight

#calcualte r.squared and plot correlation plot
wRsquared(wc$xi,wc$yi,weight=weight,order=1,plot=TRUE)
wRsquared(wc$xi,wc$yi,weight=weight,order=1,plot=TRUE,wpoints=TRUE)


#Step 2: Univariate time series fusion, using regression based prediction
xfus <- fusets(xts,yts,wc$xi,wc$yi,wc$weight)

plot2ts(xts,yts) #plot original x and y time series
plot2ts(xfus,yts) #plot fused x time series with original y time series

#################################
#B MulTiFuse - wrapper function
#################################

#Basic
xts <- dd_ndvi_cc
yts <- dd_hvhh_mt
xfus <- multifuse(xts,yts,optimize=TRUE,plot=TRUE,alpha=0.1)
xfus
plot2ts(xfus[[1]],yts)
plot2ts(xts,yts)


#set paramter
ewf = 2 #fixed ewf in case no regression weight optimization is done
optimize=TRUE #optimize regression weight
max_ewf = 2 #maximum ewf to be optimized
steps = 0.1 #optimization steps
order = 1 #regression order
plot = TRUE #plot weight optimization plot and regression plot

xfus <- multifuse(dd_ndvi_cc,dd_hvhh_mt,ewf=ewf,optimize=optimize,max_ewf=max_ewf,steps=steps,order=order,plot=TRUE)


#################################
#4) BFM
#################################

formula <- response ~ 1
order <- 1
start <- c(2007,1)
h <- 1
bfm <- bfastmonitor(xfus[[1]], formula = formula, start = start,order=order, h = h,history="ROC")
bfm <- bfastmonitor(xts, formula = formula, start = start,order=order, h = h,history="ROC")


bfm <- bfastmonitor_reich006(xts, formula = formula, start = start,order=order, h = h,history="ROC",hdyn=TRUE, history_obs_min=2)
bfm <- bfastmonitor_reich006(yts, formula = formula, start = start,order=order, h = h,history="ROC",hdyn=TRUE, history_obs_min=2)
bfm <- bfastmonitor_reich006(xfus[[1]], formula = formula,order=order, start = start,history="ROC", h = h, hdyn=TRUE, history_obs_min=2)

plot(bfm)
bfm$mag



