#LOAD DATA

library(reshape2)
library(raster)
library(ggplot2)

#CLIMATE EXTREME DATA
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Extremes/ExtremesSynched/Data/ClimateData/r1i1p1-2014-11-26/r1i1p1/")

#Load data as raster
#Warm spell duration index
wsdi.br=brick("wsdiETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc") # opens netcdf as R object

#TXx, Annual maximum value of daily maximum temperature
txx.br=brick("txxETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc") # opens netcdf as R object

#TX90p, Percentage of days when TX > 90th percentile 
tx90.br <- brick("tx90pETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")
#monthly values
#tx90.br <- brick("tx90pETCCDI_mon_NCEPREANALYSIS_historical_r1i1p1_194801-201112.nc")

#GSL, Growing season length
gsl.br <- brick("gslETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")

#TNN, Annual maximum value of daily minimum temperature
tnn.br <- brick("tnnETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")
#TNN, Monthly minimum value of daily minimum temperature
#tnn.br <- brick("tnnETCCDI_mon_NCEPREANALYSIS_historical_r1i1p1_194801-201112.nc")

#TN10p, Percentage of days when TN < 10th percentile
tn10p.br <- brick("tn10pETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")
#tn10p.br <- brick("tn10pETCCDI_mon_NCEPREANALYSIS_historical_r1i1p1_194801-201112.nc")

#CSDI, Cold spell duration index: Annual count of days with at least 6 consecutive days when TN < 10th percentile
csdi.br <- brick("csdiETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")

#DTR, Daily temperature range: Monthly mean difference between TX and TN
dtr.br <- brick("dtrETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")
#dtr.br <- brick("dtrETCCDI_mon_NCEPREANALYSIS_historical_r1i1p1_194801-201112.nc")

#set up crop to extent
extent.r= extent(-180, 180, -60, 75)   #extent(wsdi.r)
worldcrop<-crop(wrld_simpl, extent.r)

#crop to extent and average brick
for(r.k in 1:8){

  if(r.k==1) rb= wsdi.br
  if(r.k==2) rb= txx.br
  if(r.k==3) rb= tx90.br
  if(r.k==4) rb= gsl.br
  if(r.k==5) rb= tnn.br
  if(r.k==6) rb= tn10p.br
  if(r.k==7) rb= csdi.br
  if(r.k==8) rb= dtr.br
          
  rb= rotate(rb) 
  # rasterize output, give cells value of NAME(seas are NA)
  #crop
  worldcropr = rasterize(worldcrop, rb, field='NAME', fun='first')
  # mask random grid by worldcropr
  wsdi.br = mask(x=rb, mask=worldcropr)
  #mean
  rb= mean(rb)
  ##mean for 1st two layers
  #m <- mean(b[[1:2]])

  #save
  if(r.k==1) wsdi.r = rb
  if(r.k==2) txx.r = rb
  if(r.k==3) tx90.r = rb
  if(r.k==4) gsl.r = rb
  if(r.k==5) tnn.r = rb
  if(r.k==6) tn10p.r = rb
  if(r.k==7) csdi.r = rb
  if(r.k==8) dtr.r = rb
}

#---------------
#BIODIVERSITY
#IUCN: https://www.iucnredlist.org/resources/other-spatial-downloads
#mammals, birds, and amphibians
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/biodiversity/Richness_2021/")
div <- raster("Richness_2021.tif")
plot(div)

#project to lat lon
div.ll <-projectRaster(from = div, to= wsdi.r)
# mask random grid by worldcropr
div.r = mask(x=div.ll, mask=worldcropr)

#diversity R package: https://github.com/RS-eco/rasterSp
#------------------
#DISEASE

#zoonotic disease
#https://www.nature.com/articles/s41467-017-00923-8
#https://github.com/ecohealthalliance/hotspots2

#host parasite networks
#https://doi.org/10.1111/1365-2656.13666



#------------------

#compare patterns

#extract all to locate data
div.v= extract(div.ll, c(1:ncell(div.ll)))
inds= which(!is.na(div.v))
xys= xyFromCell(div.ll, inds)

#extract values
wsdi.v= extract(wsdi.r, inds)
txx.v= extract(txx.r, inds)
tx90.v= extract(tx90.r, inds)
gsl.v= extract(gsl.r, inds)
tnn.v= extract(tnn.r, inds)
tn10p.v= extract(tn10p.r, inds)
csdi.v= extract(csdi.r, inds)
dtr.v= extract(dtr.r, inds)

#diversity
div.v= extract(div.r, inds)

#combine
xy.dat= cbind(xys, wsdi.v, txx.v, tx90.v, gsl.v, tnn.v, tn10p.v, csdi.v, dtr.v, div.v)
xy.dat= as.data.frame(xy.dat)

dat.l= melt(xy.dat, id.vars= c("x","y","div.v"))

#plot relationships
ggplot(dat.l, aes(x=value, y=div.v, color=abs(y)))+geom_point()+
  facet_wrap(~variable, scales="free_x")

#maps
image(txx.r)
image(div.r)
image(wsdi.r)

#overlay maps

