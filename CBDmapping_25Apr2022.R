## PLOT LATITUDINAL GRADIENTS TOGETHER AND CHECK CORRELATIONS 

 ### PROCESS CLIMEX INDICES AND PRODUCE MAPS
#CLIMEX INDEX description: http://www.climdex.org/indices.html

#load packages 
library(ncdf4)
library(chron)
library(maptools)
library(fields)
library(raster)
library(scales) #for transparency
library(rgeos)
library(reshape2)

#source functions
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Extremes/ExtremesSynched/Analysis/ExtremeInd/")
source("ImageScale.R")

#Read data from:
#http://www.cccma.ec.gc.ca/data/climdex/climdex.shtml  ##NCEP reanalysis

###OTHER DATA SOURCES
#http://www.metoffice.gov.uk/hadobs/hadex2/index.html
#http://www.climdex.org/

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Extremes/ExtremesSynched/Data/ClimateData/r1i1p1-2014-11-26/r1i1p1/")

#===================================================
#Mask data to land
data(wrld_simpl)

#Load data as raster
#Warm spell duration index
wsdi.r=raster("wsdiETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc") # opens netcdf as R object
wsdi.r= rotate(wsdi.r) 

#crop to extent
extent.r= extent(-180, 180, -89.5, 89.5)   #extent(wsdi.r)
worldcrop<-crop(wrld_simpl, extent.r)
# rasterize output, give cells value of NAME(seas are NA)
worldcropr = rasterize(worldcrop,wsdi.r, field='NAME', fun='first')
# mask random grid by worldcropr
wsdi.r = mask(x=wsdi.r, mask=worldcropr)
#extract NAs is alignment correspondid to ncdf package
mask.r= as.matrix(wsdi.r)
nas<- is.na(mask.r)
nast= t(nas) 
nast= nast[,94:1]
#====================================================
#Load indices using ncdf package

#Warm spell duration index
wsdi.nc <- nc_open("wsdiETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc", write=TRUE) # opens netcdf as R object
# show information about the structure of the data

lat = ncvar_get(wsdi.nc, "lat")
lon = ncvar_get(wsdi.nc,"lon")
time = ncvar_get(wsdi.nc,"time") 
time_bnds= ncvar_get(wsdi.nc,"time_bnds")
wsdi= ncvar_get(wsdi.nc,"wsdiETCCDI")
yrs= 1948:2011
nc_close(wsdi.nc)

#TXx, Annual maximum value of daily maximum temperature
txx.nc <- nc_open("txxETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")
#TXx, Monthly maximum value of daily maximum temperature
#txx.nc <- nc_open("txxETCCDI_mon_NCEPREANALYSIS_historical_r1i1p1_194801-201112.nc")
print(txx.nc)
txx= ncvar_get(txx.nc,"txxETCCDI")
nc_close(txx.nc)

#TX90p, Percentage of days when TX > 90th percentile 
#tx90.nc <- nc_open("tx90pETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")
tx90.nc <- nc_open("tx90pETCCDI_mon_NCEPREANALYSIS_historical_r1i1p1_194801-201112.nc")
print(tx90.nc)
tx90= ncvar_get(tx90.nc,"tx90pETCCDI")
nc_close(tx90.nc)

#GSL, Growing season length
gsl.nc <- nc_open("gslETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")
print(gsl.nc)
gsl= ncvar_get(gsl.nc,"gslETCCDI")
nc_close(gsl.nc)

#TNN, Annual maximum value of daily minimum temperature
tnn.nc <- nc_open("tnnETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")
#TNN, Monthly minimum value of daily minimum temperature
#tnn.nc <- nc_open("tnnETCCDI_mon_NCEPREANALYSIS_historical_r1i1p1_194801-201112.nc")
print(tnn.nc)
tnn= ncvar_get(tnn.nc,"tnnETCCDI")
nc_close(tnn.nc)

#TN10p, Percentage of days when TN < 10th percentile
#tn10p.nc <- nc_open("tn10pETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")
tn10p.nc <- nc_open("tn10pETCCDI_mon_NCEPREANALYSIS_historical_r1i1p1_194801-201112.nc")
print(tn10p.nc)
tn10p= ncvar_get(tn10p.nc,"tn10pETCCDI")
nc_close(tn10p.nc)

#CSDI, Cold spell duration index: Annual count of days with at least 6 consecutive days when TN < 10th percentile
csdi.nc <- nc_open("csdiETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")
print(csdi.nc)
csdi= ncvar_get(csdi.nc,"csdiETCCDI")
nc_close(csdi.nc)

#DTR, Daily temperature range: Monthly mean difference between TX and TN
dtr.nc <- nc_open("dtrETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc")
#dtr.nc <- nc_open("dtrETCCDI_mon_NCEPREANALYSIS_historical_r1i1p1_194801-201112.nc")
print(dtr.nc)
dtr= ncvar_get(dtr.nc,"dtrETCCDI")
nc_close(dtr.nc)

#------------------
#READ TMEAN DATA
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Extremes/ExtremesSynched/Data/ClimateData/")
tmp.dat=read.csv("NCEP_Tmean.csv", header=T)
names(tmp.dat)[2:3]= c("lat","Tmean")
#restrict to lat.lad -38 to 82.85
tmp.dat= tmp.dat[which(tmp.dat$lat>-38 & tmp.dat$lat<82.85),]

#---------------------------------------

#Reorder to correspond to standard global map
lat= rev(lat)
lon[which(lon>180)]= lon[which(lon>180)]-360
ord= order(lon)
lon=lon[ord]
lat.all=lat

wsdi= wsdi[ord,length(lat):1,]
txx= txx[ord,length(lat):1,]
tx90= tx90[ord,length(lat):1,]
gsl= gsl[ord,length(lat):1,]
tnn= tnn[ord,length(lat):1,]
tn10p= tn10p[ord,length(lat):1,]
csdi= csdi[ord,length(lat):1,]
dtr= dtr[ord,length(lat):1,]

#SET start and end
start.base=(1961-1948)*12+1
end.base=(1990-1948)*12+12

start.rec=(1991-1948)*12+1
end.rec=(2010-1948)*12+12

#--------------------------

#FIX GSL>365
dims= dim(gsl)
gsl1=gsl
for(d in 1:dims[3]){
  mat1= gsl[,,d]
  mat1[which(gsl[,,d]>365)] <- 365
  gsl1[,,d]=mat1  
}
gsl=gsl1

#FIGURE 1, MAP, LATITUDE, DELTA
######### FIGURE 1

#set up plot
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")

file<-paste("CLIMEXind_Fig1.tiff" ,sep="", collapse=NULL)
#pdf(file,height = 10, width = 7)
tiff("Fig2.tiff", width = 7, height = 10, units = 'in', res = 300, compression = 'lzw')
par(cex.lab=1.2)

var.names=c("TXx (?C)", "TNn (?C)", "DTR (?C)", "GSL (days)")

layout(matrix(c(1:16),4,4,byrow=TRUE), widths=c(8,1,4,2), heights=c(1,1,1,1), FALSE)
#layout.show(16) 


for(var.k in 1:4){
  
  if(var.k==1) var=txx
  if(var.k==2) var=tnn
  if(var.k==3) var=dtr
  if(var.k==4) var=gsl #Growing season length
  
  #--------------
  #MASK TO LAND
   
  dims= dim(var)
  var.mask= var
  for(d in 1:dims[3]){
    mat1= var[,,d]
    mat1[which(nast == TRUE)] <- NA
    var.mask[,,d]=mat1  
  }
  
  #CUT OFF ANTARCTICA
  var=var.mask[,19:94,]
  lat=lat.all[19:94]
  
  #-----------------------
  #Calculate latitudinal gradient
  
 # if(var.k>2){ #yearly data
  var.base= var[,,14:43] #1961 to 1990 
  var.rec=  var[,,44:63] #1991 to 2010
#  }
  
 # if(var.k<3){ #monthly data
#    var.base= var[,,start.base:end.base] #1961 to 1990 
#    var.rec=  var[,,start.rec:end.rec] #1991 to 2010
#  }
  
  #BASE
  var= var.base
  var.ave= apply(var, MARGIN=c(1,2), FUN="mean", na.rm = TRUE) 
  var.ave[sapply(var.ave,is.nan)] = NA
  var.ave.plot=var.ave
  
  #----------------
  #ESTIMATE LAND AREA
  count= function (x) length(na.omit(x))
  land= apply(var.ave, MARGIN=2, FUN="count")
  lat.withland= which(land>10)
  
  #Restrict latitudinal gradient to land>10
  var.ave= var.ave[,lat.withland]
  var.rec= var.rec[,lat.withland,]
  lat.land= lat[lat.withland]
  
  #-------------------------
  #var.lat.mean= colMeans(var.ave, na.rm=TRUE)
  #var.lat.min= apply(var.ave, FUN="min", MARGIN=2, na.rm=TRUE)
  #var.lat.max= apply(var.ave, FUN="max", MARGIN=2, na.rm=TRUE)
  var.lat.mean= apply(var.ave, FUN="median", MARGIN=2, na.rm=TRUE)
  var.lat.min= apply(var.ave, FUN="quantile", MARGIN=2, na.rm=TRUE, probs=0.25)
  var.lat.max= apply(var.ave, FUN="quantile", MARGIN=2, na.rm=TRUE, probs=0.75)
  
  var.lat.mean[is.nan(var.lat.mean)] = NA
  var.lat.min[is.infinite(var.lat.min)] = NA
  var.lat.max[is.infinite(var.lat.max)] = NA
  
  no.na= which(!is.na(var.lat.mean))
  lat.nona= lat.land [no.na]
  var.lat.mean= na.omit(var.lat.mean)
  var.lat.min= na.omit(var.lat.min)
  var.lat.max= na.omit(var.lat.max)
  
  var.all= c(var.lat.min, var.lat.max)
  
  #-------------------------
  #RECENT
  var= var.rec
  var.ave.rec= apply(var, MARGIN=c(1,2), FUN="mean", na.rm = TRUE) 
  var.ave.rec[sapply(var.ave.rec,is.nan)] = NA
  
  #sumoflog= function(x) sum(log(na.omit(x))) 
  
  #var.lat.mean.rec= colMeans(var.ave.rec, na.rm=TRUE)
  #var.lat.min.rec= apply(var.ave.rec, FUN="min", MARGIN=2, na.rm=TRUE)
  #var.lat.max.rec= apply(var.ave.rec, FUN="max", MARGIN=2, na.rm=TRUE)
  var.lat.mean.rec= apply(var.ave.rec, FUN="median", MARGIN=2, na.rm=TRUE)
  var.lat.min.rec= apply(var.ave.rec, FUN="quantile", MARGIN=2, na.rm=TRUE, probs=0.25)
  var.lat.max.rec= apply(var.ave.rec, FUN="quantile", MARGIN=2, na.rm=TRUE, probs=0.75)
  
  var.lat.mean.rec[is.nan(var.lat.mean.rec)] = NA
  var.lat.min.rec[is.infinite(var.lat.min.rec)] = NA
  var.lat.max.rec[is.infinite(var.lat.max.rec)] = NA
  
  no.na= which(!is.na(var.lat.mean.rec))
  lat.nona.rec= lat.land[no.na]
  var.lat.mean.rec= na.omit(var.lat.mean.rec)
  var.lat.min.rec= na.omit(var.lat.min.rec)
  var.lat.max.rec= na.omit(var.lat.max.rec)
  
  var.all.rec= c(var.lat.min.rec, var.lat.max.rec)
  
  #-------------------------
  #DELTA
  var.ave.delta= var.ave.rec - var.ave
  
  var.ave.delta[sapply(var.ave.delta,is.nan)] = NA
  
  #var.lat.mean.delta= colMeans(var.ave.delta, na.rm=TRUE)
  var.lat.mean.delta= apply(var.ave.delta, FUN="median", MARGIN=2, na.rm=TRUE)
  var.lat.min.delta= apply(var.ave.delta, FUN="quantile", MARGIN=2, na.rm=TRUE, probs=0.25)
  var.lat.max.delta= apply(var.ave.delta, FUN="quantile", MARGIN=2, na.rm=TRUE, probs=0.75)
  
  var.lat.mean.delta[is.nan(var.lat.mean.delta)] = NA
  var.lat.min.delta[is.infinite(var.lat.min.delta)] = NA
  var.lat.max.delta[is.infinite(var.lat.max.delta)] = NA
  
  no.na= which(!is.na(var.lat.mean.delta))
  lat.nona.delta= lat.land[no.na]
  var.lat.mean.delta= na.omit(var.lat.mean.delta)
  var.lat.min.delta= na.omit(var.lat.min.delta)
  var.lat.max.delta= na.omit(var.lat.max.delta)
  
  var.all.delta= c(var.lat.min.delta, var.lat.max.delta)
  
  #*********************************************************
  
  #MAP BASELINE
  n=20 #number breaks
  par(mar=c(3,3,1,1), bty="o") 
  image( lon,lat, var.ave.plot, col=terrain.colors(n-1), ylab="Latitude (?)", xlab="Longitude (?)", tck=0.02, mgp=c(1, 0, 0))
  plot(wrld_simpl, add = TRUE) 
  #variable label
  #mtext(line=1, side=3, var.names[var.k], outer=F) 
  
  #Add tropics
  abline(h=23.5)
  abline(h=-23.5)
  
  #LEGEND BASELINE
  #par(oma=c( 0,0,0,0.1))# reset margin to be much smaller.
  #plot(1, type="n", axes=F, xlab="", ylab="")
  #image.plot( lon,lat, var.ave, col=terrain.colors(n-1), graphics.reset=FALSE, add=TRUE, horizontal=FALSE, legend.only=TRUE, smallplot=c(0,0.2,0.2,0.8), legend.width=1)
  par(mar=c(3,1,1,1))
  image.scale(var.ave, col = terrain.colors(n-1), horiz=FALSE, xlab="", ylab=var.names[var.k], mgp=c(1, 0.5, 0))
  
  #LAT GRAD
  par(mar=c(3,2,1,1))
  
  xlims= range(c(var.all))
  plot(var.lat.mean, lat.nona, type="l", xlim=xlims, ylim=range(lat), xlab=var.names[var.k], ylab="Latitude (0)", tck=0.02, mgp=c(1, 0, 0)) #,  yaxt='n'
  polygon(c(var.lat.min,rev(var.lat.max)),c(lat.nona,rev(lat.nona)),col="gray", border=NA)
  points(var.lat.mean, lat.nona, type="l")
  #Add recent
  #points(var.lat.mean.rec, lat.nona, type="l", col="red")
  #Add tropics
  abline(h=23.5)
  abline(h=-23.5)
  
  #ADD Tmean
  if(var.k %in% 1:2) points(tmp.dat$Tmean, tmp.dat$lat, type="l", lty="dashed")
  
  #MAP CHANGE
  #LAT GRAD
  par(mar=c(3,0.1,1,1))
  xlims= range(c(var.all.delta))
  delta.lab= bquote(Delta~.(var.names[var.k]))
  plot(var.lat.mean.delta, lat.nona.delta, type="l", ylim=range(lat), xlim=xlims, xlab=delta.lab, ylab="Latitude (0)", yaxt='n', , tck=0.02, mgp=c(1, 0, 0))
  polygon(c(var.lat.min.delta,rev(var.lat.max.delta)),c(lat.nona.delta,rev(lat.nona.delta)),col="gray", border=NA)
  points(var.lat.mean.delta, lat.nona.delta, type="l")
  #Add tropics
  abline(h=23.5)
  abline(h=-23.5)
  #zero line
  abline(v=0)
  
  #ADD Tmean
  if(var.k %in% 1:2) points(tmp.dat$tm.anom, tmp.dat$lat, type="l", lty="dashed")
  
  #STORE DATA
  if(var.k==1){
    Txx.lat.mean=var.lat.mean
    Txx.lat.min= var.lat.min
    Txx.lat.max= var.lat.max
    Txx.ave.plot=var.ave.plot
    Txx.lat.mean.delta=var.lat.mean.delta
    Txx.lat.min.delta=var.lat.min.delta
    Txx.lat.max.delta=var.lat.max.delta
    
  }
  if(var.k==2){
    Tnn.lat.mean=var.lat.mean
    Tnn.lat.min= var.lat.min
    Tnn.lat.max= var.lat.max
    Tnn.ave.plot=var.ave.plot
    Tnn.lat.mean.delta=var.lat.mean.delta
    Tnn.lat.min.delta=var.lat.min.delta
    Tnn.lat.max.delta=var.lat.max.delta
  }
  
} #end loop variables


dev.off()
#*******************************

######### FIGURE 3
#set up plot
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")
tiff("Fig3.tiff", width = 9, height = 6, units = 'in', res = 500, compression = 'lzw')

#file<-paste("Fig3.tiff" ,sep="", collapse=NULL)
#pdf(file,height = 6, width = 9)
par(cex.lab=1.2, bty="o")

var.names=c("WSI (%)", "TX90p (%)","CSI (%)","TN10p (%)", "GSL (days)","DTR (?C)")

layout(matrix(1:10,2,5,byrow=TRUE), widths=c(7,1,3,2,2,7,1,3,2,2), heights=c(1,1,1,1,1,1,1,1,1,1), FALSE)
#layout.show(10)  

for(var.k in 1:4){
  
  if(var.k==1) var=wsdi
  if(var.k==2) var=tx90
  if(var.k==3) var=csdi
  if(var.k==4) var= tn10p
  
  #--------------
  #MASK TO LAND
  
  dims= dim(var)
  var.mask= var
  for(d in 1:dims[3]){
    mat1= var[,,d]
    mat1[which(nast == TRUE)] <- NA
    var.mask[,,d]=mat1  
  }
  
  #CUT OFF ANTARCTICA
  var=var.mask[,19:94,]
  lat=lat.all[19:94]
  
  #-----------------------
  #Calculate latitudinal gradient
  
  var.base= var[,,14:43] #1961 to 1990 
  var.rec=  var[,,44:63] #1991 to 2010
  
  #BASE
  var= var.base
  
  #vary metrics across years 
  ppos= function(x) sum(na.omit(x)>0)/length(na.omit(x)) 
  
  var.ave= apply(var, MARGIN=c(1,2), FUN="mean", na.rm = TRUE) 
  
  if(var.k %in% c(1,3))
  {  
    var.ave= apply(var, MARGIN=c(1,2), FUN="ppos")*100 
  }
  
  var.ave[sapply(var.ave,is.nan)] = NA
  var.ave[sapply(var.ave,is.infinite)] = NA
  var.ave.plot=var.ave
  
  #sumoflog= function(x) sum(log(na.omit(x))) 
  
  #----------------
  #ESTIMATE LAND AREA
  count= function (x) length(na.omit(x))
  land= apply(var.ave, MARGIN=2, FUN="count")
  lat.withland= which(land>10)
  
  #Restrict latitudinal gradient to land>10
  var.ave= var.ave[,lat.withland]
  var.rec= var.rec[,lat.withland,]
  lat.land= lat[lat.withland]
  #-----------------
  
  #var.lat.mean= colMeans(var.ave, na.rm=TRUE)
  #var.lat.min= apply(var.ave, FUN="min", MARGIN=2, na.rm=TRUE)
  #var.lat.max= apply(var.ave, FUN="max", MARGIN=2, na.rm=TRUE)
  var.lat.mean= apply(var.ave, FUN="median", MARGIN=2, na.rm=TRUE)
  var.lat.min= apply(var.ave, FUN="quantile", MARGIN=2, na.rm=TRUE, probs=0.25)
  var.lat.max= apply(var.ave, FUN="quantile", MARGIN=2, na.rm=TRUE, probs=0.75)
  
  var.lat.mean[is.nan(var.lat.mean)] = NA
  var.lat.min[is.infinite(var.lat.min)] = NA
  var.lat.max[is.infinite(var.lat.max)] = NA
  
  no.na= which(!is.na(var.lat.mean))
  lat.nona= lat.land [no.na]
  var.lat.mean= na.omit(var.lat.mean)
  var.lat.min= na.omit(var.lat.min)
  var.lat.max= na.omit(var.lat.max)
  
  var.all= c(var.lat.min, var.lat.max)
  
  #-------------------------
  #RECENT
  var= var.rec
  
  var.ave.rec= apply(var, MARGIN=c(1,2), FUN="mean", na.rm = TRUE) 
  
  if(var.k %in% c(1,3))
  {  
    # var.ave.rec= apply(var, MARGIN=c(1,2), FUN="mean", na.rm = TRUE) 
    #  var.ave.rec= apply(var, MARGIN=c(1,2), FUN="max", na.rm = TRUE) 
    var.ave.rec= apply(var, MARGIN=c(1,2), FUN="ppos") *100
  }
  
  var.ave.rec[sapply(var.ave.rec,is.nan)] = NA
  var.ave.rec[sapply(var.ave.rec,is.infinite)] = NA
  
  #-------------------------
  #DELTA
  var.ave.delta= var.ave.rec - var.ave
  
  var.ave.delta[sapply(var.ave.delta,is.nan)] = NA
  
  #var.lat.mean.delta= colMeans(var.ave.delta, na.rm=TRUE)
  var.lat.mean.delta= apply(var.ave.delta, FUN="median", MARGIN=2, na.rm=TRUE)
  var.lat.min.delta= apply(var.ave.delta, FUN="quantile", MARGIN=2, na.rm=TRUE, probs=0.25)
  var.lat.max.delta= apply(var.ave.delta, FUN="quantile", MARGIN=2, na.rm=TRUE, probs=0.75)
  
  var.lat.mean.delta[is.nan(var.lat.mean.delta)] = NA
  var.lat.min.delta[is.infinite(var.lat.min.delta)] = NA
  var.lat.max.delta[is.infinite(var.lat.max.delta)] = NA
  
  no.na= which(!is.na(var.lat.mean.delta))
  lat.nona.delta= lat.land[no.na]
  var.lat.mean.delta= na.omit(var.lat.mean.delta)
  var.lat.min.delta= na.omit(var.lat.min.delta)
  var.lat.max.delta= na.omit(var.lat.max.delta)
  
  var.all.delta= c(var.lat.min.delta, var.lat.max.delta)
  
  #*********************************************************
  if(var.k %in% c(1,3,5)){
    
    
    #MAP BASELINE
    n=20 #number breaks
    par(mar=c(4,4,1,1), mgp=c(2, 1, 0))
    image( lon,lat, var.ave.plot, col=terrain.colors(n-1), ylab="Latitude (?)", xlab="Longitude (?)", tck=0.02, mgp=c(1, 0, 0))
    plot(wrld_simpl, add = TRUE) 
    #variable label
    #mtext(line=1, side=3, var.names[var.k], outer=F) 
    
    #Add tropics
    abline(h=23.5)
    abline(h=-23.5)
    
    #LEGEND BASELINE
    #par(oma=c( 0,0,0,0.1))# reset margin to be much smaller.
    #plot(1, type="n", axes=F, xlab="", ylab="")
    #image.plot( lon,lat, var.ave, col=terrain.colors(n-1), graphics.reset=FALSE, add=TRUE, horizontal=FALSE, legend.only=TRUE, smallplot=c(0,0.2,0.2,0.8), legend.width=1)
    par(mar=c(3,1,1,1), mgp=c(2, 1, 0))
    image.scale(var.ave, col = terrain.colors(n-1), horiz=FALSE, xlab="", ylab="")
    
    #LAT GRAD
    par(mar=c(4,2,1,1),mgp=c(2, 1, 0))
    
    xlims= range(c(var.all))
    plot(var.lat.mean, lat.nona, type="l", xlim=xlims, ylim=range(lat), xlab=var.names[var.k], ylab="Latitude (?)", tck=0.02, mgp=c(1, 0, 0) ) #,  yaxt='n'
    polygon(c(var.lat.min,rev(var.lat.max)),c(lat.nona,rev(lat.nona)),col="gray", border=NA)
    points(var.lat.mean, lat.nona, type="l")
    #Add recent
    #points(var.lat.mean.rec, lat.nona, type="l", col="red")
    #Add tropics
    abline(h=23.5)
    abline(h=-23.5)
    
    #MAP CHANGE
    #LAT GRAD
    par(mar=c(4,0.5,1,1), mgp=c(2, 1, 0))
    xlims= range(c(var.all.delta))
    delta.lab= bquote(Delta~.(var.names[var.k]))
    plot(var.lat.mean.delta, lat.nona.delta, type="l", ylim=range(lat), xlim=xlims, xlab=delta.lab, ylab="Latitude (?)", yaxt='n', tck=0.02, mgp=c(1, 0, 0))
    polygon(c(var.lat.min.delta,rev(var.lat.max.delta)),c(lat.nona.delta,rev(lat.nona.delta)),col="gray", border=NA)
    points(var.lat.mean.delta, lat.nona.delta, type="l")
    #Add tropics
    abline(h=23.5)
    abline(h=-23.5)
    #zero line
    abline(v=0)
  } # END CHECK VAR.K FOR LEFT COLUMN
  
  #-------------------------------------------
  if(var.k %in% c(2,4,6)){
    #LAT GRAD
    par(mar=c(4,2,1,1))
    
    #MAP CHANGE
    #LAT GRAD
    par(mar=c(4,0.5,1,1))
    xlims= range(c(var.all.delta))
    delta.lab= bquote(Delta~.(var.names[var.k]))
    plot(var.lat.mean.delta, lat.nona.delta, type="l", ylim=range(lat), xlim=xlims, xlab=delta.lab, ylab="Latitude (?)", yaxt='n', tck=0.02, mgp=c(1, 0, 0))
    polygon(c(var.lat.min.delta,rev(var.lat.max.delta)),c(lat.nona.delta,rev(lat.nona.delta)),col="gray", border=NA)
    points(var.lat.mean.delta, lat.nona.delta, type="l")
    #Add tropics
    abline(h=23.5)
    abline(h=-23.5)
    #zero line
    abline(v=0)
  } # END CHECK VAR.K FOR RIGHT COLUMN
  
  
} #end loop variables


dev.off()

#======================================
