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
rndvi80 <-dropLayer(rndvi80,c(1:37))
rndvi70 <-dropLayer(rndvi70,c(1:37))
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

###Getting single pixel
cell <- 2901
#x[2]<-8056630
hv_ts <- bfastts(as.vector(rhv[cell]),as.Date(getZ(rhv)),type=c("irregular"))
ndvi_ts <- bfastts(as.vector(rndvi90[cell]),as.Date(getZ(rndvi)),type=c("irregular"))
ndvi_ts[193] <- 0.8210577
ndvi_ts[1296] <- NA
plot2ts(ndvi_ts,hv_ts,lab_ts1="Landsat NDVI",lab_ts2="ALOS PALSAR HV [dB]")

save(ndvi_ts,hv_ts,rndvi,rndvi90,rhv,rhh,rhvhh,stableforest,loggedforest,file="fiji.Rdata")

################

