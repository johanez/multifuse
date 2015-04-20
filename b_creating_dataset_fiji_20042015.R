require(raster)
#require(zoo)
#require(bfast)  ## this loads the bfast package
#require(tseries)
#require(PreFAST)
#require(caret)
#require(psychometric)
#require(bfastSpatial)
#require(bayts)

#load all pixel based functions
source("Functions.24012014.R")
source("Functions.correlation.24012014.R")

##################
#load raster data
load("/media//DATA2/reich006/Scripts/R/sar_landsat_ts_fusion_3/data/Data_75_72_2000_2012__750_1400_1000_1650.RData")
#load("data/Data_75_72_2000_2012__750_1400_1000_1650_north.RData")
e_new <- c(563995,572495,8054965,8060465)
e_new <- c(565000,567000,8056000,8058000)
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
hv_dates

#restrict to observation period (2005)
rndvi <- dropLayer(rndvi,c(1:37))
#rndvi <- dropLayer(rndvi,c(81:96))
rndvi90 <- dropLayer(rndvi90,c(1:37))
#rndvi90 <- dropLayer(rndvi90,c(81:96))
ndvi_dates <- ndvi_dates[38:133]
#ndvi_dates <- ndvi_dates[1:80]

############
rhh <- setZ(rhh_mt,hh_dates,name='dates')
rhv <- setZ(rhv_mt,hv_dates,name='dates')
rhvhh <- setZ(rhvhh_mt,hv_dates,name='dates')
rndvi <- setZ(rndvi,ndvi_dates,name='dates')
rndvi90 <- setZ(rndvi90,ndvi_dates,name='dates')



#Reference data
plot(ylog)

#Refernece data
plot(ylog)
plot(yplant)
plot(ync)
ylogP <- ylog+((yperiod/10)*2.5)-0.25
loggedforest <- ylogP  
stableforest <- ync
plot(stableforest)
plot(loggedforest)

#layover shadow map
rhv[ls_map==1]<-NA
rhh[ls_map==1]<-NA
rhvhh[ls_map==1]<-NA
rndvi[ls_map==1]<-NA
rndvi90[ls_map==1]<-NA
stableforest[ls_map==1]<-NA
loggedforest[ls_map==1]<-NA

save(rndvi,rndvi90,rhv,rhh,rhvhh,stableforest,loggedforest,file="fiji.Rdata")

##############

#check percentage Missing data in NDVI time series
plot(calc.brick.percNA(rndvi,cores=4))

#Increase missing data pecentage
rndvi95 <- get.brick.tsMD_random(rndvi,0.9,cores=4)
names(rndvi95) <- names(rndvi)
plot(calc.brick.percNA(rndvi95,cores=4))

#Reference data
plot(ylog)



################

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
