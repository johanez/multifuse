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
#2) Weighted correlation
#################################
#wc <- ts.correlation_A3.plot(dd_ndvi, dd_hvhh_mt,maxdt= NULL,itype="two-ways",interpolation_x="linear",interpolation_y="linear",order_x=1,order_y=1,maxbreaks_x=5,maxbreaks_y=5,name_x="NDVI",name_y="HVHH")
#wc <- ts.correlation_A3(dd_ndvi_cc, dd_hvhh_mt,maxdt=NULL,itype="two-ways",interpolation_x="linear",interpolation_y="linear",order_x=1,order_y=1,maxbreaks_x=5,maxbreaks_y=5)
wc <- regressionWeight(dd_ndvi_cc,dd_hvhh_mt)

#plot original overlapping ts
plot2ts(wc$x,wc$y)
#plot interpolated overlatpping ts
plot2ts(wc$xi,wc$yi)


# 3 Calculate magnitude weights
#plot
#calc.Rsquare_weight.plot(wc$xi,wc$yi,weight=wc$w,order=1)
#calculate magnitude using exponential weight function (ewf)
#ewf=0 (no weight)
#ewf=1 (simple weight)
#ewf=2 (squared weight)
ewf=1
wmago <- wc$w^ewf
wmago <- wmago * 1/sum(wmago)
wmago
#calc.Rsquare_weight.plot(wc$xi,wc$yi,weight=wmago,order=1)
wRsquared(wc$xi,wc$yi,weight=wmago,order=1,plot=TRUE)

#calculate correlation measures
vx <- as.vector(wc$xi)
vy <- as.vector(wc$yi)
fit=lm(na.omit(vy)~na.omit(vx),weight=wc$w^2)
summary(fit)
fit$pvalue

####################
#Weight optimization
opt <- opt_wt.plot(wc$xi,wc$yi,wc$w,max_wt_exp=10,steps=0.1,order=1)
opt <- optimizeWeight(wc$xi,wc$yi,wc$w,max_ewf=10,steps=0.1,order=1,plot=TRUE)
opt

#################################
#3) Fusion
#################################

ewf <- 2
#if weight optimization 
max_ewf <- 5
wt_opt <- TRUE

dd_ndvi_ccm <- ts.fusion.correlation_A3.plot(dd_ndvi_cc,dd_hvhh_mt,wt_type = "magnitude",wt_opt=wt_opt,wt_exp=ewf,wt_opt_max_exp=max_ewf,itype="two-ways",name_x="NDVI",name_y="HVHH")
#fix(ts.fusion.correlation_A3.plot)

plot.2ts(dd_ndvi_cc,dd_hvhh_mt,name_x="NDVI",name_y="HVHH (mt)",points=TRUE)
plot.2ts(dd_ndvi_ccm,dd_hvhh_mt,name_x="Fused NDVI",name_y="HVHH (mt)",points=TRUE)


#################################
#4) BFM
#################################

formula <- response ~ 1
order <- 1
start <- c(2008,1)
h <- 0.25

bfm <- bfastmonitor_reich006(dd_ndvi_cc, formula = formula, start = start,order=order, h = h,history="ROC",hdyn=TRUE, history_obs_min=2)
bfm <- bfastmonitor_reich006(dd_ndvi_ccm, formula = formula, start = start,order=order, h = h,history="ROC",hdyn=TRUE, history_obs_min=2)
bfm <- bfastmonitor_reich006(dd_hvhh_mt, formula = formula,order=order, start = start,history="ROC", h = h, hdyn=TRUE, history_obs_min=2)

plot(bfm)




