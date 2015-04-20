#require(raster)
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
#load fiji raster data
load("fiji.Rdata")

#show reference data
plot(loggedforest)
plot(stableforest,legend=FALSE,add=TRUE)

#extract pixel time series
plot(rhv,3)
cell <- click(rhv, n=1, cell=TRUE)[,1]

#cell<-2901
#create time series using bfastts
hv <- bfastts(as.vector(rhv[cell]),as.Date(getZ(rhv)),type=c("irregular"))
ndvi <- bfastts(as.vector(rndvi[cell]),as.Date(getZ(rndvi)),type=c("irregular"))
ndvi90 <- bfastts(as.vector(rndvi90[cell]),as.Date(getZ(rndvi90)),type=c("irregular"))

#1) Plot SAR data over NDVI
plot2ts(ndvi,hv,lab_ts1="Landsat NDVI [MD=org]",lab_ts2="PALSAR HV [dB]")
plot2ts(ndvi90,hv,lab_ts1="Landsat NDVI [MD=95]",lab_ts2="PALSAR HV [dB]")

################################
#A MulTiFuse - single steps
#################################
xts <- ndvi90
yts <- hv

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
xts <- ndvi90
yts <- hv
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



