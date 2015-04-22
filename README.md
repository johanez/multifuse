# multifuse

Package to apply the Multi-sensor Time series Fusion (MulTiFuse) approach to two univariat time series (e.g. optical and SAR time series). 
The MulTiFuse approach has been published in [Reiche et al. 2015](http://www.sciencedirect.com/science/article/pii/S0034425714003885), 
where it was used to fuse univariate Landsat NDVI and ALOS PALSAR L-band SAR backscatter time series.

### The MulTiFuse approach
A detailed description of the MulTiFuse approach is provided in [Reiche et al. 2015](http://www.sciencedirect.com/science/article/pii/S0034425714003885). A brief description to fuse two univariate remote sensing time series (e.g. SAR and optical time series) is given below:


First, a weighted time series correlation is performed. To maximise the statistical significance of the correlation and to take exceptional cases into account, correlation weight optimization is done before the relationship of the two time series is modelled through a weighted regression analysis. 
The optimized regressionmodel is utilized in a second step for regression-based prediction of time series observation to fuse the two time series.


![fig2_20052014](https://cloud.githubusercontent.com/assets/6399980/7251311/77775dc4-e82a-11e4-8b6b-083cc9051fb8.jpg)
Figure 1. Schematic overview of MulTiFuse to fuse a univariate optical and SAR time series [Reiche et al. 2015](http://www.sciencedirect.com/science/article/pii/S0034425714003885).


### References

Reiche, J., Verbesselt, J., Hoekman, D. H. & Herold, M. (2015): Fusing Landsat and SAR time series to detect deforestation in the tropics. Remote Sensing of Environment. 156, 276-293. DOI: 10.1016/j.rse.2014.10.001. http://www.sciencedirect.com/science/article/pii/S0034425714003885 

# Using the multifuse package

### Installation
The package can be installed directly from github using devtools
```r
library(devtools)
install_github('jreiche/multifuse')
```
### Apply MulTiFuse
````r
## load multifuse package
library(multifuse)

## load example data
data(tsexample)
## plot example time series
plot2ts(ndvi,hv,lab_ts1="Landsat NDVI",lab_ts2="ALOS PALSAR HV [dB]")

## apply multifuse
fus <- multifuse(ndvi,hv,optimize=TRUE,ewf_max=2,ewf_steps=0.1, plot=TRUE)
## get fused time series from multifuse output
ndvi_fused<- fus[[1]]

##plot original ndvi time series and fused time series
plot2ts(ndvi,ndvi_fused,lab_ts1="Landsat NDVI",lab_ts2="Fused Landsat NDVI",ylim_ts1=c(0.3,0.9),ylim_ts2=c(0.3,0.9))

````

### Example data
(A) (Single pixel data) Landsat NDVI (2005 - 2012) and ALOS PALSAR HV (2007 - 2010) example time series
```r
## load example data
data(tsexample)

## plot example time series
plot2ts(ndvi,hv,lab_ts1="Landsat NDVI",lab_ts2="ALOS PALSAR HV [dB]")
```
(B) (Raster data) Landsat NDVI (2005 - 2012), ALOS PALSAR HV, HH and HVHH-ratio (2007 - 2010) time series for a pinus caribea planatation in Fiji. Landsat NDVI data is provided with original per pixel missing data and with 90 percent per pixel missing data. Three-monthly harvesting reference data are provided. For a detailed data description refer to [Reiche et al. 2015](http://www.sciencedirect.com/science/article/pii/S0034425714003885). 
```r
## load example data
data(fiji)

## show reference data
plot(loggedforest)
plot(stableforest,legend=FALSE,add=TRUE)

## extract pixel time series
#option 1: define cell  
cell<-2901
#option 2: select cell  
plot(rhv,3)
cell <- click(rhv, n=1, cell=TRUE)[,1]

#create time series using bfastts (bfast package)
hv <- bfastts(as.vector(rhv[cell]),as.Date(getZ(rhv)),type=c("irregular"))
hh <- bfastts(as.vector(rhh[cell]),as.Date(getZ(rhh)),type=c("irregular"))
hvhh <- bfastts(as.vector(rhvhh[cell]),as.Date(getZ(rhvhh)),type=c("irregular"))
ndvi <- bfastts(as.vector(rndvi[cell]),as.Date(getZ(rndvi)),type=c("irregular"))
ndvi90 <- bfastts(as.vector(rndvi90[cell]),as.Date(getZ(rndvi90)),type=c("irregular"))

#plot time series
plot2ts(ndvi,hv,lab_ts1="Landsat NDVI [MD=org]",lab_ts2="PALSAR HV [dB]")
plot2ts(ndvi90,hv,lab_ts1="Landsat NDVI [MD=90]",lab_ts2="PALSAR HV [dB]")
```

