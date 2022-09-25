#LOAD DATA

library(reshape2)
library(raster)
library(ggplot2)
library(viridisLite)

library(maptools)
data(wrld_simpl)

hdir=getwd()

#set up crop to extent
extent.r= extent(-180, 180, -60, 75)   #extent(wsdi.r)
rgeos::set_RGEOS_CheckValidity(2L)
worldcrop<-crop(wrld_simpl, extent.r)

#CLIMATE EXTREME DATA
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Extremes/ExtremesSynched/Data/ClimateData/r1i1p1-2014-11-26/r1i1p1/")

#Load data as raster
#Warm spell duration index
wsdi.br=brick("wsdiETCCDI_yr_NCEPREANALYSIS_historical_r1i1p1_1948-2011.nc") # opens netcdf as R object

#crop
rb= wsdi.br
rb= rotate(rb) 
worldcropr = rasterize(worldcrop, rb, field='NAME', fun='first')
wsdi.r = rb

#-------------------
#Climate Hotspots
#https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021EF002027
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/climate/MiaoHotspot/")
chot <- raster("SED_ssp245_2080-2099.nc")
chot= flip(t(chot), 1)
chot=(flip(chot, 2))

crs(chot) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
res(chot)=1
#project to lat lon
chot.ll <-projectRaster(from = chot, to= wsdi.r)
## mask random grid by worldcropr
#chot.r = mask(x=chot.ll, mask=worldcropr)

chot.ll<-crop(chot.ll, extent.r)

#---------------
#BIODIVERSITY

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
chot.v= extract(chot.ll, c(1:ncell(chot.ll)))
inds= which(!is.na(chot.v))
chot.v= chot.v[inds]

xys= xyFromCell(chot.ll, inds)

#extract values
bii.v= extract(bii.r, inds)
hpar.v= extract(hpar.r, inds)

#combine
xy.dat= cbind(xys, chot.v, hpar.v, bii.v) 
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

clim.cor <- cor(xy.dat[,c(3:5)])
corrplot(clim.cor, method = "ellipse", type = "lower")

xy.plot= xy.dat[,c("chot.v", "bii.v","hpar.v") ]

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")
pdf("CBDcorrelations.pdf",height = 8, width = 8)
ggpairs(xy.plot) 
dev.off()

#----------------- 
#plot relationships
dat.l= melt(xy.dat, id.vars= c("x","y","div.v"))

#plot relationships
ggplot(dat.l, aes(x=value, y=div.v, color=abs(y)))+geom_point()+
  facet_wrap(~variable, scales="free_x")

#3 way
xy.dat= as.data.frame(cbind(xys, chot.v, bii.v, hpar.v)) #ext.v
xy.dat$bii.threat= 1-xy.dat$bii.v
  
ggplot(xy.dat, aes(x=bii.v, y=hpar.v, color=chot.v))+geom_point()

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")
pdf("CBDscatter.pdf",height = 8, width = 8)
ggplot(xy.dat, aes(x=chot.v, y=hpar.v, color=bii.v))+geom_point(alpha=0.5)+
  theme_classic(base_size = 20)+
  ylab("host-parasite interactions")+xlab("climate hotspot")+
  scale_color_viridis_c("Biodiversity risk (1-bii)")
dev.off()

#overlay maps
#scale to max and add
cbd.int= chot/cellStats(chot, stat='max') +(1-bii.r/cellStats(bii.r, stat='max')) +
  hpar.r/cellStats(hpar.r, stat='max')

bii.risk= 1-bii.r

#=========================

threat.stack= stack(txx.r, 1-bii.r, hpar.r) 
plotRGB(threat.stack)

#3way color map
#http://matthewkling.github.io/media/colormap/

#devtools::install_github("matthewkling/colormap")
library(colormap)
library(fields)

#3 way
xy.dat= as.data.frame(cbind(xys, chot.v, bii.v, hpar.v)) 
xy.dat$bii.threat= 1-xy.dat$bii.v
xy.dat= na.omit(xy.dat)

xy.dat$cbd.int= xy.dat$chot.v/max(xy.dat$chot.v)+xy.dat$bii.threat/max(xy.dat$bii.threat)+xy.dat$hpar.v/max(xy.dat$hpar.v)

# map color to the climate variables
xy.dat$colors <- colors3d(xy.dat[,c(3,4,5)])

plot(xy.dat$x, xy.dat$y, col=xy.dat$colors)

risk.cm=ggplot(xy.dat, aes(x= x, y=y, fill=colors))+geom_tile()+
  scale_fill_identity()

quilt.plot(xy.dat$x, xy.dat$y, xy.dat$cbd.int, 
           add.legend=TRUE, col=viridis(n=20))

#write out
setwd(hdir)
write.csv(xy.dat,"threatmapdata.csv")

#3D
#https://www.rayshader.com/

#install.packages("rayshader")
library(rayshader)

tmap= ggplot(xy.dat, aes(x= x, y=y, fill=bii.threat))+geom_tile()+
  scale_fill_viridis_c()

plot_gg(tmap)

#-----
#plot together
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")
pdf("CBDmaps.pdf",height = 10, width = 6)
par(mfrow=c(4,1), mar=c(2,2,1,0.5) )

#maps
plot(chot, main="climate standard euclidean distance", xlim=c(-150,160), ylim=c(-60,80))
plot(bii.risk, main="biodiversity risk (1-bii)", xlim=c(-150,160), ylim=c(-60,80))
plot(hpar.r, main="predicted host-parasite interactions", xlim=c(-150,160), ylim=c(-60,80))
#plot(cbd.int, main="CBD overlap (sum with each each scaled to 1)", xlim=c(-150,160), ylim=c(-60,80))
dev.off()


