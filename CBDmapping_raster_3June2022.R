#LOAD DATA

library(reshape2)
library(raster)
library(ggplot2)
library(viridisLite)

library(maptools)
data(wrld_simpl)

#CLIMATE CHANGE
#http://climexp.knmi.nl/plot_atlas_form.py
#mean rcp45 temperature 2081-2100 minus 1986-2005 Jan-Dec AR5 CMIP5 subset
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/climate/")
clim.change <- raster("diff_tas_Amon_onemean_rcp45_000_2081-2100_minus_1986-2005_mon1_ave12_withsd.nc")

#global mean temperatures
#download from https://data.ceda.ac.uk/badc/cru/data/cru_ts/cru_ts_4.02/data/tmp
tmp <- raster("cru_ts4.02.1901.2017.tmp.dat.nc")

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
rgeos::set_RGEOS_CheckValidity(2L)
worldcrop<-crop(wrld_simpl, extent.r)

#crop
rb= wsdi.br
rb= rotate(rb) 
worldcropr = rasterize(worldcrop, rb, field='NAME', fun='first')

#crop to extent
# mask random grid by worldcropr
clim.change= rotate(clim.change)
clim.change.r=clim.change  #mask(x=clim.change, mask=worldcropr)

tmp.r=tmp  #mask(x=tmp, mask=worldcropr)

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

#Garcia- dimensions of climate change
#DOI: 10.1126/science.1247579
#change in probability of local climate extremes
#future data 2069-2099
#exceed 95th percentile

#Change in TX90p, Percentage of days when TX > 90th percentile 
rb= rotate(tx90.br) 
# rasterize output, give cells value of NAME(seas are NA)
#crop
worldcropr = rasterize(worldcrop, rb, field='NAME', fun='first')
# mask random grid by worldcropr
rb = mask(x=rb, mask=worldcropr)
#mean
tx90.init= mean(rb[[1:49]])
tx90.rec= mean(rb[[50:64]])
#change in percent
tx90.dif= (tx90.rec-tx90.init)/tx90.init
plot(tx90.dif)

#probability of at least one record breaking extreme per year
#https://www.nature.com/articles/s41558-021-01092-9
#https://data.iac.ethz.ch/Fischer_et_al_2021_RecordExtremes/figures/

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/climate/")
ext <- raster("fig3b.nc")
ext= rotate(ext)

#project to lat lon
ext.ll <-projectRaster(from = ext, to= wsdi.r)
# mask random grid by worldcropr
ext.r = mask(x=ext.ll, mask=worldcropr)

#---------------
#BIODIVERSITY
#IUCN: https://www.iucnredlist.org/resources/other-spatial-downloads
#mammals, birds, and amphibians
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/biodiversity/Richness_2021/")
div <- raster("Richness_2021.tif")

#project to lat lon
div.ll <-projectRaster(from = div, to= wsdi.r)
# mask random grid by worldcropr
div.r = mask(x=div.ll, mask=worldcropr)

# #diversity R package: https://github.com/RS-eco/rasterSp
# 
# #BII
# #https://data.nhm.ac.uk/dataset/global-map-of-the-biodiversity-intactness-index-from-newbold-et-al-2016-science
# setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/biodiversity/BII/")
# bii <- raster("lbii.asc")
# 
# crs(bii) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
# 
# #project to lat lon
# bii.ll <-projectRaster(from = bii, to= wsdi.r)
# # mask random grid by worldcropr
# bii.r = mask(x=bii.ll, mask=worldcropr)
#save projected
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/out/")
#saveRDS(bii.r, "bii.rds")
bii.r= readRDS("bii.rds")
#------------------
#DISEASE

#zoonotic disease
#https://www.nature.com/articles/s41467-017-00923-8
#https://github.com/ecohealthalliance/hotspots2

#host parasite networks
#https://doi.org/10.1111/1365-2656.13666

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/disease/HostParasite_raster/")
hpar= readRDS("phylogenyModel_raster.rds")
#also "affinityModel_raster.rds","observed_links_raster.rds","phylogenyModel_raster.rds", "combinedModel_raster.rds"
#observed
#hpar.obs= readRDS("observed_links_raster.rds")

#project to lat lon
hpar.ll <-projectRaster(from = hpar, to= wsdi.r)
# mask random grid by worldcropr
hpar.r = mask(x=hpar.ll, mask=worldcropr)
  
#------------------

#compare patterns

#extract all to locate data
div.v= extract(div.ll, c(1:ncell(div.ll)))
inds= which(!is.na(div.v))
xys= xyFromCell(div.ll, inds)

#extract values
clim.change.v= extract(clim.change.r, inds)
tmp.v= extract(tmp.r, inds)

wsdi.v= extract(wsdi.r, inds)
txx.v= extract(txx.r, inds)
tx90.v= extract(tx90.r, inds)
gsl.v= extract(gsl.r, inds)
tnn.v= extract(tnn.r, inds)
tn10p.v= extract(tn10p.r, inds)
csdi.v= extract(csdi.r, inds)
dtr.v= extract(dtr.r, inds)

#extremes
ext.v= extract(ext.r, inds)

tx90dif.v= extract(tx90.dif, inds)

#diversity
div.v= extract(div.r, inds)
bii.v= extract(bii.r, inds)

#disease
hpar.v= extract(hpar.r, inds)

#combine
xy.dat= cbind(xys, clim.change.v, tmp.v, wsdi.v, txx.v, tx90.v, gsl.v, tnn.v, tn10p.v, csdi.v, dtr.v, div.v, ext.v, tx90dif.v,hpar.v, bii.v) 
xy.dat= as.data.frame(xy.dat)

#write out
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/out/")
#write.csv(xy.dat,"pts.csv")
#xy.dat= read.csv("pts.csv")

#-----------------
#CORRELATIONS
library(corrplot)
#library(Cairo)
#library(psych)
library(GGally) #can add groups https://r-charts.com/correlation/ggpairs/

clim.cor <- cor(xy.dat[,c(3:12,15)])
corrplot(clim.cor, method = "ellipse", type = "lower")

#try out climate pcs
xy.omit= na.omit(xy.dat[,c(5:7,9:12)])
pc=princomp(xy.omit, cor=FALSE, fix_sign=FALSE)
pc$loadings

#xy.plot= cbind(xy.dat[,c("x","y", "div.v", "bii.v","hpar.v") ] ,pc=pc$scores[,1])
#first pc is just txx.v and tnn.v that are highly correlated 
xy.plot= xy.dat[,c("txx.v", "tnn.v", "div.v", "bii.v","hpar.v") ]

ggpairs(xy.plot) 

#----------------- 
#plot relationships
dat.l= melt(xy.dat, id.vars= c("x","y","div.v"))

#plot relationships
ggplot(dat.l, aes(x=value, y=div.v, color=abs(y)))+geom_point()+
  facet_wrap(~variable, scales="free_x")

#3 way
xy.dat= as.data.frame(cbind(xys, bii.v, ext.v, hpar.v))
xy.dat$bii.threat= 1-xy.dat$bii.v
  
ggplot(xy.dat, aes(x=bii.v, y=hpar.v, color=ext.v))+geom_point()

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")
pdf("CBDscatter.pdf",height = 8, width = 8)
ggplot(xy.dat, aes(x=ext.v, y=hpar.v, color=bii.threat))+geom_point(alpha=0.5)+
  theme_classic(base_size = 20)+
  ylab("host-parasite interactions")+xlab("probability of annual heat extreme")+
  scale_color_viridis_c("Biodiversity risk (1-bii)")
dev.off()

library(rgl)
plot3d(x=ext.v, y=1-xy.dat$bii.v, z=hpar.v)

#overlay maps
#scale to max and add
#cbd.int= ext.r/cellStats(ext.r, stat='max') +(1-bii.r/cellStats(bii.r, stat='max')) +
#  hpar.r/cellStats(hpar.r, stat='max')

#use max
cbd.int= txx.r/cellStats(txx.r, stat='max') +(1-bii.r/cellStats(bii.r, stat='max')) +
  hpar.r/cellStats(hpar.r, stat='max')

bii.risk= 1-bii.r

#plot together
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")
pdf("CBDmaps.pdf",height = 10, width = 6)
par(mfrow=c(4,1), mar=c(2,2,1,0.5) )

#maps
#plot(tx90.dif)
plot(ext.r, main="probability of annual record breaking heat extreme", xlim=c(-150,160), ylim=c(-60,80))
plot(bii.risk, main="biodiversity risk (1-bii)", xlim=c(-150,160), ylim=c(-60,80))
plot(hpar.r, main="predicted host-parasite interactions", xlim=c(-150,160), ylim=c(-60,80))
plot(cbd.int, main="CBD overlap (sum with each each scaled to 1)", xlim=c(-150,160), ylim=c(-60,80))

dev.off()
