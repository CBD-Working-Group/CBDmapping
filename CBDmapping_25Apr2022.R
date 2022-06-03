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

#---------------------------------------

#Reorder to correspond to standard global map
lat= rev(lat)
lon[which(lon>180)]= lon[which(lon>180)]-360
ord= order(lon)
lon=lon[ord]

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

######### FIGURE 2

#set up plot
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Extremes\\Figs\\GCBFigs\\")
#setwd("C:\\Users\\Buckley\\Google Drive\\BuckleySynch\\Extremes\\Figs\\IndFigs1\\")
#setwd("C:\\Users\\lbuckley\\Google Drive\\BuckleySynch\\Extremes\\Figs\\IndFigs\\")

tiff("Fig4.tiff", width = 6, height = 6, units = 'in', res = 500, compression = 'lzw')
#file<-paste("CLIMEXind_Fig2.pdf" ,sep="", collapse=NULL)
#cairo_pdf(file,height = 6, width = 6)

layout(matrix(1:6,2,3,byrow=TRUE), widths=c(5,4,4), heights=c(1,1,1), FALSE)
par(mar=c(3,4,1,1), bty="o", tck=0.02, mgp=c(1, 0, 0))

xlims= range(-60,60)

#INSECT PLOT
plot(Txx.lat.mean, lat.nona, type="l", xlim=xlims, ylim=range(lat), xlab="Temperature (?C)", ylab="Latitude (?)") #,  yaxt='n'
polygon(c(Txx.lat.min,rev(Txx.lat.max)),c(lat.nona,rev(lat.nona)),col="gray", border=NA)
points(Txx.lat.mean, lat.nona, type="l")

polygon(c(Tnn.lat.min,rev(Tnn.lat.max)),c(lat.nona,rev(lat.nona)),col="gray", border=NA)
points(Tnn.lat.mean, lat.nona, type="l")

#add mean
points(tmp.dat$Tmean, tmp.dat$lat, type="l", col="black", lty="dashed")

#Add tropics
abline(h=23.5)
abline(h=-23.5)

#ADD THERMAL TOLERANCE DATA
#points(hoff.up$MetricValue, hoff.up$Latitude, pch=-9658, col=alpha("black", 0.5) )
#points(hoff.low$MetricValue, hoff.low$Latitude, pch=-9668, col=alpha("black", 0.5) ) 

points(hoff.CTmax$MetricValue, hoff.CTmax$Latitude, pch=24, bg=alpha("black", 0.5), col=alpha("black", 0.5) )
points(hoff.CTmin$MetricValue, hoff.CTmin$Latitude, pch=25, bg=alpha("black", 0.5), col=alpha("black", 0.5) ) 

points(hoff.ULT$MetricValue, hoff.ULT$Latitude, pch=2, col=alpha("black", 0.5) )
points(hoff.LLT$MetricValue, hoff.LLT$Latitude, pch=6, col=alpha("black", 0.5) ) 

hoff.low= rbind(hoff.LLT[,c("Latitude", "MetricValue")],  hoff.CTmin[,c("Latitude", "MetricValue")])
hoff.up= rbind(hoff.ULT[,c("Latitude", "MetricValue")],  hoff.CTmax[,c("Latitude", "MetricValue")])
hoff.low.n= subset(hoff.low, hoff.low$Latitude>=0)
hoff.low.s= subset(hoff.low, hoff.low$Latitude<0)
hoff.up.n= subset(hoff.up, hoff.up$Latitude>=0)
hoff.up.s= subset(hoff.up, hoff.up$Latitude<0)

for(g in 1:4){

if(g==1) dat= hoff.low.n[order(hoff.low.n$Latitude),]
if(g==2) dat= hoff.low.s[order(hoff.low.s$Latitude),]
if(g==3) dat= hoff.up.n[order(hoff.up.n$Latitude),]
if(g==4) dat= hoff.up.s[order(hoff.up.s$Latitude),]

mod1= lm(dat$"MetricValue" ~dat$Latitude)
pred1= predict(mod1)

segments(pred1[1], dat[1,1], pred1[length(pred1)], dat[nrow(dat),1], lwd=1,lty="dotted", col="red")
}


#-------------------
#OTHER ECTO PLOT
par(mar=c(3,0.5,1,1))
plot(Txx.lat.mean, lat.nona, type="l", xlim=xlims, ylim=range(lat), xlab="Temperature (?C)", ylab="Latitude (?)",  yaxt='n') #
polygon(c(Txx.lat.min,rev(Txx.lat.max)),c(lat.nona,rev(lat.nona)),col="gray", border=NA)
points(Txx.lat.mean, lat.nona, type="l")

polygon(c(Tnn.lat.min,rev(Tnn.lat.max)),c(lat.nona,rev(lat.nona)),col="gray", border=NA)
points(Tnn.lat.mean, lat.nona, type="l")

#add mean
points(tmp.dat$Tmean, tmp.dat$lat, type="l", col="black", lty="dashed")

#Add tropics
abline(h=23.5)
abline(h=-23.5)

#ADD THERMAL TOLERANCE DATA
#divide reptiles and amphibians
amp.d= sun.dat.lt[sun.dat.lt$class=="Amphibia",]
rep.d= sun.dat.lt[sun.dat.lt$class=="Reptilia",]
for(m in 1:nrow(rep.d)) points(c(rep.d$tmin[m],rep.d$tmax[m]), rep(rep.d$lat[m],2), type="l", col=alpha("blue", 0.5),lty="dotted")
for(m in 1:nrow(amp.d)) points(c(amp.d$tmin[m],amp.d$tmax[m]), rep(amp.d$lat[m],2), type="l", col=alpha("green", 0.5),lty="dotted")

#divide reptiles and amphibians
amp.d= sun.dat[sun.dat$class=="Amphibia",]
rep.d= sun.dat[sun.dat$class=="Reptilia",]
#plot as lines
for(m in 1:nrow(rep.d)) points(c(rep.d$tmin[m],rep.d$tmax[m]), rep(rep.d$lat[m],2), type="l", col=alpha("blue", 0.5))
for(m in 1:nrow(amp.d)) points(c(amp.d$tmin[m],amp.d$tmax[m]), rep(amp.d$lat[m],2), type="l", col=alpha("green", 0.5))

#----------------

amp.n= subset(amp.d, amp.d$lat>=0)
amp.s= subset(amp.d, amp.d$lat<0)
rep.n= subset(rep.d, rep.d$lat>=0)
rep.s= subset(rep.d, rep.d$lat<0)
cols=c("green","green","blue","blue")

for(g in c(1,3,4)){
  
  if(g==1) dat= amp.n[order(amp.n$lat),]
  if(g==2) dat= amp.s[order(amp.s$lat),]
  if(g==3) dat= rep.n[order(rep.n$lat),]
  if(g==4) dat= rep.s[order(rep.s$lat),]
  
  mod1= lm(dat$tmin ~dat$lat)
  pred1= predict(mod1)
  segments(pred1[1], dat[1,"lat"], pred1[length(pred1)], dat[nrow(dat),"lat"], lwd=1, col=cols[g],lty="dotted")

  mod1= lm(dat$tmax ~dat$lat)
  pred1= predict(mod1)
  segments(pred1[1], dat[1,"lat"], pred1[length(pred1)], dat[nrow(dat),"lat"], lwd=1, col=cols[g],lty="dotted")
}

#------------------------
#ANIMAL DIVERSITY 
mam.dat$rel.div= mam.dat$div / max(mam.dat$div)
amp.dat$rel.div= amp.dat$AmpSR / max(amp.dat$AmpSR)
turt.dat$rel.div= turt.dat$Nspecies / max(turt.dat$Nspecies)
ant.dat$rel.div= ant.dat$Total.number.of.genera / max(ant.dat$Total.number.of.genera)
bird.div$rel.div= bird.div$div / max(bird.div$div)
rep.div$rel.div= rep.div$ReptileSR / max(rep.div$ReptileSR)

par(mar=c(3,0.5,1,1))
plot(mam.dat$rel.div, mam.dat$lats,  type="l", col="purple", lty="solid", xlim= range(0, 1),ylim=range(lat), ylab="Latitude (?)", xlab="Diversity", lwd=1, yaxt='n')

points(rep.div$rel.div, rep.div$lat, type="l", col="blue", lwd=1)
points(turt.dat$rel.div, turt.dat$latitude, type="l", col="blue", lwd=1, lty="dashed")
points(ant.dat$rel.div, ant.dat$Latitude, type="l", col="black", lwd=1)

points(amp.dat$rel.div, amp.dat$Lat, type="l", col="green", lwd=1)
points(bird.div$rel.div, bird.div$lats, type="l", col="orange", lwd=1)

legend("topright",legend=c("Mammals", "Birds","Amphibians", "Reptiles",  "Turtles", "Ant genera"),col=c("purple","orange","green", "blue","blue", "black"),lty=c("solid","solid","solid", "solid","dashed", "solid"), bty="n")

#Add tropics
abline(h=23.5)
abline(h=-23.5)

#----------------------
#ENDO PLOT  
#BIRDS
par(mar=c(3,4,1,1))

plot(Txx.lat.mean, lat.nona, type="l", xlim=xlims, ylim=range(lat), xlab="Temperature (?C)", ylab="Latitude (?)")
polygon(c(Txx.lat.min,rev(Txx.lat.max)),c(lat.nona,rev(lat.nona)),col="gray", border=NA)
points(Txx.lat.mean, lat.nona, type="l")

polygon(c(Tnn.lat.min,rev(Tnn.lat.max)),c(lat.nona,rev(lat.nona)),col="gray", border=NA)
points(Tnn.lat.mean, lat.nona, type="l")

#add mean
points(tmp.dat$Tmean, tmp.dat$lat, type="l", col="black", lty="dashed")

#Add tropics
abline(h=23.5)
abline(h=-23.5)

#ADD THERMAL TOLERANCE DATA
#points(mammal.dat$LTNZ...C., mammal.dat$Latitude...N., pch=1, col=alpha("purple", 0.5))
#points(mammal.dat$UTNZ...C., mammal.dat$Latitude...N., pch=1, col=alpha("purple", 0.5))
#points(bird.dat$LTNZ...C., bird.dat$Latitude...N., pch=1, col=alpha("green", 0.5))
#points(bird.dat$UTNZ...C., bird.dat$Latitude...N., pch=1, col=alpha("green", 0.5))

#plot as lines
for(m in 1:nrow(bird.dat)) points(c(bird.dat$LTNZ...C.[m],bird.dat$UTNZ...C.[m]), rep(bird.dat$Latitude...N.[m],2), type="l", col=alpha("orange", 0.5))


#----------------
bird.n= subset(bird.dat, bird.dat$Latitude...N.>=0)
bird.s= subset(bird.dat, bird.dat$Latitude...N.<0)

for(g in 1:2){
  
  if(g==1) dat= bird.n[order(bird.n$Latitude...N.),]
  if(g==2) dat= bird.s[order(bird.s$Latitude...N.),]
  
  mod1= lm(dat$UTNZ...C. ~dat$Latitude...N.)
  pred1= predict(mod1)
  segments(pred1[1], dat[1,"Latitude...N."], pred1[length(pred1)], dat[nrow(dat),"Latitude...N."], lwd=1, lty="dotted", col="orange")
  
  mod1= lm(dat$LTNZ...C. ~dat$Latitude...N.)
  pred1= predict(mod1)
  segments(pred1[1], dat[1,"Latitude...N."], pred1[length(pred1)], dat[nrow(dat),"Latitude...N."], lwd=1, lty="dotted", col="orange")
}

#---------------
#MAMMALS
par(mar=c(3,0.5,1,1))

plot(Txx.lat.mean, lat.nona, type="l", xlim=xlims, ylim=range(lat), xlab="Temperature (?C)", ylab="Latitude (?)",  yaxt='n')
polygon(c(Txx.lat.min,rev(Txx.lat.max)),c(lat.nona,rev(lat.nona)),col="gray", border=NA)
points(Txx.lat.mean, lat.nona, type="l")

polygon(c(Tnn.lat.min,rev(Tnn.lat.max)),c(lat.nona,rev(lat.nona)),col="gray", border=NA)
points(Tnn.lat.mean, lat.nona, type="l")

#add mean
points(tmp.dat$Tmean, tmp.dat$lat, type="l", col="black", lty="dashed")

#Add tropics
abline(h=23.5)
abline(h=-23.5)

#plot as lines
for(m in 1:nrow(mammal.dat)) points(c(mammal.dat$LTNZ...C.[m],mammal.dat$UTNZ...C.[m]), rep(mammal.dat$Latitude...N.[m],2), type="l", col=alpha("purple", 0.5))

#----------------
mammal.n= subset(mammal.dat, mammal.dat$Latitude...N.>=0)
mammal.s= subset(mammal.dat, mammal.dat$Latitude...N.<0)

for(g in 1:2){
  
  if(g==1) dat= mammal.n[order(mammal.n$Latitude...N.),]
  if(g==2) dat= mammal.s[order(mammal.s$Latitude...N.),]
  
  mod1= lm(dat$UTNZ...C. ~dat$Latitude...N.)
  pred1= predict(mod1)
  segments(pred1[1], dat[1,"Latitude...N."], pred1[length(pred1)], dat[nrow(dat),"Latitude...N."], lwd=1,lty="dotted", col="purple")
  
  mod1= lm(dat$LTNZ...C. ~dat$Latitude...N.)
  pred1= predict(mod1)
  segments(pred1[1], dat[1,"Latitude...N."], pred1[length(pred1)], dat[nrow(dat),"Latitude...N."], lwd=1,lty="dotted", col="purple")
}

#--------------------------
#crops and population
par(mar=c(3,0.5,1,1))

#aggregate
crop.area$lat.bin= round(crop.area$lats)
crop.area= aggregate(crop.area, list(crop.area$lat.bin), FUN="mean")

crop.yield$lat.bin= round(crop.yield$lats)
crop.yield= aggregate(crop.yield, list(crop.yield$lat.bin), FUN="mean")

pop$lat.bin= round(pop$lats)
pop= aggregate(pop, list(pop$lat.bin), FUN="mean")

crop.area$rel.area= crop.area$lat.grad / max(crop.area$lat.grad, na.rm=TRUE)
crop.yield$rel.yield= crop.yield$crop.lat.grad / max(crop.yield$crop.lat.grad)
pop$rel.yield= pop$pop.lat.grad / max(pop$pop.lat.grad)

plot(crop.area$rel.area,crop.area$lats,  type="l",lwd=1,lty="dashed", col="black", xlim= range(0, 1),ylim=range(lat), ylab="Latitude (?)", xlab="Metric", yaxt='n') #proportion
points( crop.yield$rel.yield, crop.yield$lats, type="l", col="black")
points( pop$rel.yield,pop$lats, type="l", col="grey")

legend("topright",legend=c("Population", "% Crop Area", "Crop Yield"),lty=c("solid","dashed","solid"),col=c("grey", "black","black"), bty="n")


#Add tropics
abline(h=23.5)
abline(h=-23.5)

dev.off()

######### FIGURE 3
#set up plot
#setwd("C:\\Users\\Buckley\\Google Drive\\BuckleySynch\\Extremes\\Figs\\IndFigs1\\")
#setwd("C:\\Users\\lbuckley\\Google Drive\\BuckleySynch\\Extremes\\Figs\\IndFigs\\")
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

#Text statistics for WSI

ws= cbind(var.lat.mean, lat.nona)
colMeans(ws[which(lat.nona>23),])

ws=cbind(var.lat.mean.delta, lat.nona.delta)

#------------------------------------------
#FIGURE 5
#CORRELATIONS
library(corrplot)
library(Cairo)
library(psych)

cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
#---------------

#set up plot
#setwd("C:\\Users\\Buckley\\Google Drive\\BuckleySynch\\Extremes\\Figs\\IndFigs1\\")
#setwd("C:\\Users\\lbuckley\\Google Drive\\BuckleySynch\\Extremes\\Figs\\IndFigs\\")

file<-paste("Fig5.tiff" ,sep="", collapse=NULL)
#cairo_pdf(file,height = 6, width = 10)
tiff(file, width = 10, height = 6, units = 'in', res = 500, compression = 'lzw')

layout(matrix(1:2,1,2,byrow=TRUE), widths=c(0.8,1.0), heights=c(1,1), FALSE)
#layout.show(2)

vars= list(txx,tnn,wsdi,tx90,csdi, tn10p,gsl,dtr)

for(var.k in 1:8){
  var= vars[[var.k]]
  
    var.base= var[,,14:43] #1961 to 1990 
    var.rec=  var[,,44:63] #1991 to 2010
  
  var.ave= apply(var[,,14:43], MARGIN=c(1,2), FUN="mean", na.rm = TRUE) 
  var.ave.rec= apply(var[,,44:63], MARGIN=c(1,2), FUN="mean", na.rm = TRUE) 
  
  var.ave[sapply(var.ave,is.nan)] = NA
  var.ave.rec[sapply(var.ave.rec,is.nan)] = NA
  var.delta= var.ave.rec - var.ave
  
  if(var.k==1) {delta.mat= c(var.delta); var.mat=c(var.ave)}
  if(var.k>1) {delta.mat= cbind(delta.mat,c(var.delta)); var.mat= cbind(var.mat,c(var.ave))}
  
} #end loop var.k

#assign names
colnames(delta.mat)= c("  TXx","  TNn","  WSI ","  TX90p"," CSI", "  TN10p","  GSL","  DTR")
colnames(var.mat)= c("  TXx","  TNn","  WSI ","  TX90p"," CSI", "  TN10p","  GSL","  DTR")

var.mat.ord=  var.mat[,c("  TXx","  TNn","  GSL","  DTR","  WSI ","  TX90p"," CSI", "  TN10p")]
bg.cols= c(rep("darkgray",10), rep("lightgray",4), rep("tan",1), rep("lightgray",4), rep("tan",2), rep("lightgray",4), rep("tan",3), rep("lightgray",4), rep("tan",4))

M <- cor(var.mat.ord[,c(1:5,7)])
res1 <- cor.mtest(var.mat.ord[,c(1:5,7)], 0.95)
cp=corrplot(M, method = "ellipse",p.mat = res1[[1]], sig.level = 0.05, order = "original", type = "upper", tl.pos="lt", bg=bg.cols, tl.col = "black")
#write.csv(M, "BaselineCor.csv")

bg.cols= c(rep("darkgray",4), rep("lightgray",2),rep("darkgray",3), rep("lightgray",2),rep("darkgray",2), rep("lightgray",2),rep("darkgray",1), rep("lightgray",2), rep("tan",10) )

delta.mat.ord= delta.mat[,dimnames(cp)[[1]] ]
M <- cor(delta.mat.ord)
res1 <- cor.mtest(delta.mat.ord, 0.95)
corrplot(M, method = "ellipse", p.mat = res1[[1]], sig.level = 0.05, order = "original", type = "lower", add=TRUE, tl.pos="n", bg=bg.cols, tl.col = "black", cl.pos="n")
#write.csv(M, "DeltaCor.csv")

#----------------------------------
## DELTAS
delta.mat.ord=  delta.mat[,c("  TXx","  TNn","  GSL","  DTR","  WSI "," CSI","  TX90p", "  TN10p")]
colnames(delta.mat.ord)=c("\u0394 TXx" ,"\u0394 TNn","\u0394 GSL","\u0394 DTR","\u0394 WSI","\u0394 CSI","\u0394 TX90p","\u0394 TN10p")
 
var.mat.ord=  var.mat[,c("  TXx","  TNn","  GSL","  DTR","  WSI "," CSI")]

bg.cols= c(rep("darkgray",10), rep("lightgray",4), rep("tan",1), rep("lightgray",4), rep("tan",2), rep("lightgray",4), rep("tan",2), rep("lightgray",4), rep("tan",2))

M<-cor(var.mat.ord, delta.mat.ord)
M<-corr.test(var.mat.ord, delta.mat.ord, adjust='none')
corrplot(M$r, method = "ellipse", p.mat = M$p, sig.level = 0.05, order = "original", type = "upper", tl.pos="lt", bg= bg.cols, tl.col = "black")
#write.csv(M, "BaselineDeltaCor.csv")
#cairo_pdf(pdf.file,family="Arial Unicode MS")

dev.off()

corr.test(var.mat.ord, delta.mat.ord, adjust='none')
#------------------------------------------
#FIGURE S1
#THERMODYNAMIC

#Estimate thermodynamic mean
tmp.dat$Ta.therm= exp(-E/(k*(tmp.dat$tm.rec+273.15)))
tmp.dat$TaAnom.therm= exp(-E/(k*(tmp.dat$tm.anom+273.15)))

#set up plot
setwd("C:\\Users\\Buckley\\Google Drive\\BuckleySynch\\Extremes\\Figs\\IndFigs1\\")
#setwd("C:\\Users\\lbuckley\\Google Drive\\BuckleySynch\\Extremes\\Figs\\IndFigs\\")

file<-paste("CLIMEXind_FigS1.pdf" ,sep="", collapse=NULL)
pdf(file,height = 5, width = 6)

var.names=c("TXx (?C)", "TNn (?C)", "DTR (?C)", "GSL (days)")

layout(matrix(c(1:8),2,4,byrow=TRUE), widths=c(6,1,4,2), heights=c(1,1,1,1), FALSE)
#layout.show(16) 


for(var.k in 1:2){
  
  if(var.k==1) var=txx.therm
  if(var.k==2) var=tnn.therm
  
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
  
  #-----------------------
  #Calculate latitudinal gradient
  
  #if(var.k>2){ #yearly data
    var.base= var[,,14:43] #1961 to 1990 
    var.rec=  var[,,44:63] #1991 to 2010
  #}
  
  #BASE
  var= var.base
  var.ave= apply(var, MARGIN=c(1,2), FUN="mean", na.rm = TRUE) 
  var.ave[sapply(var.ave,is.nan)] = NA
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
  
 
  
  #*********************************************************
  
  #MAP BASELINE
  n=20 #number breaks
  par(mar=c(3,3,1,1))
  image( lon,lat, var.ave.plot, col=terrain.colors(n-1), ylab="Latitude (?)", xlab="Longitude (?)", mgp=c(2, 1, 0))
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
  plot(var.lat.mean, lat.nona, type="l", xlim=xlims, ylim=range(lat), xlab=var.names[var.k], ylab="Latitude (0)", mgp=c(2, 1, 0)) #,  yaxt='n'
  polygon(c(var.lat.min,rev(var.lat.max)),c(lat.nona,rev(lat.nona)),col="gray", border=NA)
  points(var.lat.mean, lat.nona, type="l")
  #Add recent
  #points(var.lat.mean.rec, lat.nona, type="l", col="red")
  #Add tropics
  abline(h=23.5)
  abline(h=-23.5)
  
  #ADD Tmean
  if(var.k %in% 1:2) points(tmp.dat$Ta.therm, tmp.dat$lat, type="l", lty="dashed")
  
  #MAP CHANGE
  #LAT GRAD
  par(mar=c(3,0.1,1,1))
  xlims= range(c(var.all.delta))
  delta.lab= bquote(Delta~.(var.names[var.k]))
  plot(var.lat.mean.delta, lat.nona.delta, type="l", ylim=range(lat), xlim=xlims, xlab=delta.lab, ylab="Latitude (0)", yaxt='n', mgp=c(2, 1, 0))
  polygon(c(var.lat.min.delta,rev(var.lat.max.delta)),c(lat.nona.delta,rev(lat.nona.delta)),col="gray", border=NA)
  points(var.lat.mean.delta, lat.nona.delta, type="l")
  #Add tropics
  abline(h=23.5)
  abline(h=-23.5)
  #zero line
  abline(v=0)
  
  #ADD Tmean
  if(var.k %in% 1:2) points(tmp.dat$TaAnom.therm, tmp.dat$lat, type="l", lty="dashed")
  
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

#---------------------------
#OTHER PLOTS
### PCA
plot(var.mat)
pr= prcomp(var.mat, center = TRUE, scale. = TRUE)
summary(pr)
print(pr)

biplot(pr,xlabs=rep("o",nrow(var.mat)))

####################################################
# PLOT LATITUDINAL GRADIENTS TOGETHER

vars= list(txx,tnn,wsdi,tx90,csdi, tn10p,gsl,dtr)
lat.all=lat

for(var.k in 1:8){
  var= vars[[var.k]]
  
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
  #------------------
  
  var.ave= apply(var[,,14:43], MARGIN=c(1,2), FUN="mean", na.rm = TRUE) 
  var.ave.rec= apply(var[,,44:63], MARGIN=c(1,2), FUN="mean", na.rm = TRUE) 
  
  #----------------
  #ESTIMATE LAND AREA
  count= function (x) length(na.omit(x))
  land= apply(var.ave, MARGIN=2, FUN="count")
  lat.withland= which(land>10)
  
  #Restrict latitudinal gradient to land>10
  var.ave= var.ave[,lat.withland]
  var.ave.rec= var.ave.rec[,lat.withland]
  lat.land= lat[lat.withland]
  #-------------------
  
  #lat gradients
  #vary metrics across years 
  ppos= function(x) sum(na.omit(x)>0)/length(na.omit(x)) 
  
  var.ave1= apply(var.ave, MARGIN=c(2), FUN="mean", na.rm = TRUE) 
  var.rec1= apply(var.ave.rec, MARGIN=c(2), FUN="mean", na.rm = TRUE) 
  
  if(var.k %in% c(3,5))
  {  
    var.ave1= apply(var.ave, MARGIN=c(2), FUN="ppos") 
    var.rec1= apply(var.ave.rec, MARGIN=c(2), FUN="ppos")
  }
  
  var.ave= var.ave1
  var.rec= var.rec1
  
  var.ave[sapply(var.ave,is.nan)] = NA
  var.ave[sapply(var.ave,is.infinite)] = NA
  var.rec[sapply(var.rec,is.nan)] = NA
  var.rec[sapply(var.rec,is.infinite)] = NA
  
  var.delta= var.rec - var.ave
  #-----------------
  
  if(var.k==1) { lat.mat= var.ave; lat.mat.rec=var.rec;  lat.mat.delta=var.delta}
  if(var.k>1) {lat.mat= cbind(lat.mat,c(var.ave)); lat.mat.rec= cbind(lat.mat.rec,c(var.rec)); lat.mat.delta= cbind(lat.mat.delta,c(var.delta))}
  
} #end loop var.k

#assign names
colnames(lat.mat)= c("TXx","TNn","wsdi","TX90p","csdi", "TN10p","gsl","dtr")
colnames(lat.mat.rec)= c("TXx","TNn","wsdi","TX90p","csdi", "TN10p","gsl","dtr")
colnames(lat.mat.delta)= c("TXx","TNn","wsdi","TX90p","csdi", "TN10p","gsl","dtr")

#reorder
lat.mat= lat.mat[,c("TXx","TNn","dtr", "gsl", "wsdi","csdi", "TX90p","TN10p")]
lat.mat.rec= lat.mat.rec[,c("TXx","TNn","dtr", "gsl", "wsdi","csdi", "TX90p","TN10p")]
lat.mat.delta= lat.mat.delta[,c("TXx","TNn","dtr", "gsl", "wsdi","csdi", "TX90p","TN10p")]

#Normalize columns to 1
for(c in 1:8){
  lat.mat[,c]= (lat.mat[,c]-min(lat.mat[,c]))/(max(lat.mat[,c])- min(lat.mat[,c]))
  lat.mat.rec[,c]= (lat.mat.rec[,c]-min(lat.mat.rec[,c]))/(max(lat.mat.rec[,c])- min(lat.mat.rec[,c]))
  lat.mat.delta[,c]= (lat.mat.delta[,c])/max(abs(lat.mat.delta[,c]))
}

#plot
cols= c("red","blue","orange","green", "red","blue","salmon", "cyan3")
ltys= c("solid","solid","solid","solid", "dashed","dashed","dashed", "dashed")

par(mfrow=c(1,3),mar=c(3,3,2,1), las=1, oma=c(2,2,0,4.5), xpd=FALSE)
plot(lat.mat[,1],lat.land,  type="l", lwd=2, xlim= c(0,1), col= cols[1],lty=ltys[1], main="baseline", xlab="")
for(x in 1:6){ points( lat.mat[,x],lat.land, type="l", lwd=2, col= cols[x],lty=ltys[x] )}

plot(lat.mat.rec[,1],lat.land,  type="l", lwd=2, xlim= c(0,1), col= cols[1],lty=ltys[1], main="recent", xlab="" )
for(x in 2:8){ points( lat.mat.rec[,x],lat.land, type="l", lwd=2, col= cols[x],lty=ltys[x])}

plot(lat.mat.delta[,1],lat.land,  type="l", lwd=2, xlim= c(-1,1), col= cols[1],lty=ltys[1], main="delta", xlab="" )
for(x in 2:8){ points( lat.mat.delta[,x],lat.land, type="l", lwd=2, col= cols[x],lty=ltys[x] )}
abline(v=0, lwd=1)

par(xpd=NA)
legend(1,80,col=cols, colnames(lat.mat), lty=ltys, bty="n", lwd=2)

mtext("Latitude",line=0, side=2, outer=TRUE, las=0, cex=1.2) 
mtext("Normalized index",line=0, side=1, outer=TRUE, cex=1.2)

####################################################
#library("psych")
data(sat.act)
corr.test(sat.act)
corr.test(M, adjust = "none")
