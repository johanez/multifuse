# multifuse

Package to apply the Multi-sensor Time series Fusion (MulTiFuse) approach to two univariat time series. 
The MulTiFuse approach has been published in [Reiche et al. 2015](http://www.sciencedirect.com/science/article/pii/S0034425714003885), 
and was used to fuse univariate Landsat NDVI and ALOS PALSAR backscatter time series.

### The MulTiFuse approach
A detailed description of the MulTiFuse approach is provided in [Reiche et al. 2015](http://www.sciencedirect.com/science/article/pii/S0034425714003885). A brief description is given below:


First, a weighted timeseries correlation is performed. To maximise the statistical significance of the correlation and to take exceptional cases into account, correlation weight optimization is done before the relationship of the two time series is modelled through a weighted regression analysis. 
The optimized regressionmodel is utilized in a second step for regression-based prediction of time series observation to fuse the two time series.


![fig2_20052014](https://cloud.githubusercontent.com/assets/6399980/7251311/77775dc4-e82a-11e4-8b6b-083cc9051fb8.jpg)
Schematic overview of MulTiFuse to fuse a optical and SAR time series [Reiche et al. 2015](http://www.sciencedirect.com/science/article/pii/S0034425714003885).


### Installation
The package can be installed directly from github using devtools
```
library(devtools)
install_github('jreiche/multifuse')
```
### References

Reiche, J., Verbesselt, J., Hoekman, D. H. & Herold, M. (2015): Fusing Landsat and SAR time series to detect deforestation in the tropics. Remote Sensing of Environment. 156, 276-293. DOI: 10.1016/j.rse.2014.10.001. http://www.sciencedirect.com/science/article/pii/S0034425714003885 

# Using the multifuse package

### Loading example data
(A) Landsat NDVI (2005 - 2012) and ALOS PALSAR HV (2007 - 2010) example time series
```
## load example data
load("tsexample.rda")
xts <- ndvi #Landsat NDVI example time series
yts <- hv   #ALSO PALSAR HV example time series
```
(B) (Raster data) Landsat NDVI (2005 - 2012), ALOS PALSAR HV, HH and HVHH-ratio (2007 - 2010) time series.


## plot original time series
plot2ts(xts,yts,lab_ts1="Landsat NDVI",lab_ts2="ALOS PALSAR HV [dB]")
 
## apply multifuse
xfus <- multifuse(xts,yts,optimize=TRUE,ewf_max=2,ewf_steps=0.1, plot=TRUE)

##plot fused time series
plot2ts(xfus[[1]],yts,lab_ts1="Fused Landsat NDVI",lab_ts2="ALOS PALSAR HV [dB]")
```
