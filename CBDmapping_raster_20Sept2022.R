#LOAD DATA

library(reshape2)
library(raster)
library(ggplot2)
library(viridisLite)
library(patchwork)

library(maptools)
data(wrld_simpl)

hdir=getwd()

#set up crop to extent
extent.r= extent(-180, 180, -60, 75)   #extent(wsdi.r)
rgeos::set_RGEOS_CheckValidity(2L)
worldcrop<-crop(wrld_simpl, extent.r)

#------------------
#DISEASE

#zoonotic disease from Hahn et al.
hr <- raster("Total host richness")
hr.ll<-crop(hr, extent.r) 

#------------------

#Climate Hotspots
#https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021EF002027
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/climate/MiaoHotspot/")
chot <- raster("SED_ssp245_2080-2099.nc")
chot= flip(t(chot), 1)
chot=(flip(chot, 2))

crs(chot) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
res(chot)=1
chot.ll<-crop(chot, extent.r)

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
# bii.ll <-projectRaster(from = bii, to= hr.ll)
# # mask random grid by worldcropr
# bii.r = mask(x=bii.ll, mask=worldcropr)
#save projected
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/out/")
#saveRDS(bii.r, "bii.rds")
bii.r= readRDS("bii.rds")

#https://figshare.com/articles/dataset/Global_maps_of_Biodiversity_Intactness_Index_Sanchez-Ortiz_et_al_2019_-_bioRxiv_/7951415/1
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/biodiversity/BII2019/abundance_richness_BII_maps/")
bii2<- raster("final-rich-bii-isl-main.tif")
bii.r = crop(bii2, extent.r)

#------------------
#compare patterns

#extract all to locate data
chot.v= extract(chot.ll, c(1:ncell(chot.ll)))
inds= which(!is.na(chot.v))
chot.v= chot.v[inds]

xys= xyFromCell(chot.ll, inds)

#extract values
bii.v= extract(bii.r, xys)
hr.v= extract(hr.ll, xys)

#combine
xy.dat= cbind(xys, chot.v, hr.v, bii.v) 
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

xy.plot= xy.dat[,c("chot.v", "bii.v","hr.v") ]

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")
pdf("CBDcorrelations.pdf",height = 8, width = 8)
ggpairs(xy.plot) 
dev.off()

#----------------- 
#plot relationships
dat.l= melt(xy.dat, id.vars= c("x","y"))

#plot relationships
ggplot(dat.l, aes(x=value, y=div.v, color=abs(y)))+geom_point()+
  facet_wrap(~variable, scales="free_x")

#3 way
xy.dat= as.data.frame(cbind(xys, chot.v, bii.v, hr.v)) #ext.v
xy.dat$bii.threat= 1-xy.dat$bii.v
  
ggplot(xy.dat, aes(x=bii.v, y=hr.v, color=chot.v))+geom_point()

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")
pdf("CBDscatter.pdf",height = 8, width = 8)
ggplot(xy.dat, aes(x=chot.v, y=hr.v, color=bii.v))+geom_point(alpha=0.5)+
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

threat.stack= stack(txx.r, 1-bii.r, hr.r) 
plotRGB(threat.stack)

#3way color map
#http://matthewkling.github.io/media/colormap/

#devtools::install_github("matthewkling/colormap")
library(colormap)
library(fields)

#3 way
xy.dat= as.data.frame(cbind(xys, chot.v, bii.v, hr.v)) 
xy.dat$bii.threat= 1-xy.dat$bii.v
xy.dat= na.omit(xy.dat)

xy.dat$cbd.int= xy.dat$chot.v/max(xy.dat$chot.v)+xy.dat$bii.threat/max(xy.dat$bii.threat)+xy.dat$hr.v/max(xy.dat$hr.v)

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

# #maps
# plot(chot, main="climate standard euclidean distance", xlim=c(-150,160), ylim=c(-60,80))
# plot(bii.risk, main="biodiversity risk (1-bii)", xlim=c(-150,160), ylim=c(-60,80))
# plot(hpar.r, main="predicted host-parasite interactions", xlim=c(-150,160), ylim=c(-60,80))
# print(risk.cm)
# #plot(cbd.int, main="CBD overlap (sum with each each scaled to 1)", xlim=c(-150,160), ylim=c(-60,80))

#climate
chot.df <- as.data.frame(chot, xy=TRUE) #Convert raster to data.frame
names(chot.df)[3] <- 'climate' 
head(chot.df)

clim.plot= ggplot(data = chot.df)+
  geom_raster(mapping=aes(x=x, y=y, fill=climate))+
  scale_fill_gradientn(colours= rev(terrain.colors(10)), name='climate change', na.value="white")+
  theme_void()

#biodiversity
bii.df <- as.data.frame(bii.risk, xy=TRUE) #Convert raster to data.frame
names(bii.df)[3] <- 'biodiversity' 
head(chot.df)

bii.plot= ggplot(data = bii.df)+
  geom_raster(mapping=aes(x=x, y=y, fill=biodiversity))+
  scale_fill_gradientn(colours= rev(terrain.colors(10)), name='biodiversity risk', na.value="white")+
  theme_void()

#disease
hpar.df <- as.data.frame(hr.ll, xy=TRUE) #Convert raster to data.frame
names(hpar.df)[3] <- 'disease' 
head(hpar.df)

disease.plot= ggplot(data = hpar.df)+
  geom_raster(mapping=aes(x=x, y=y, fill=disease))+
  scale_fill_gradientn(colours= rev(terrain.colors(10)), name='host richness', na.value="white")+ #'predicted host-parasite interactions'
  theme_void()

#-----------
#Risk map
# map color to the climate variables
xy.dat$colors <- colors3d(xy.dat[,c(3,4,5)])
#names(xy.dat)[3:5]=c("climate change","biodiversity risk","host richness")

#plot(xy.dat$x, xy.dat$y, col=xy.dat$colors)

risk.cm=ggplot(xy.dat, aes(x= x, y=y, fill=colors))+geom_tile()+
  scale_fill_identity()+theme_void()

#plot legend
library(rgl)
with(xy.dat, plot3d(chot.v, bii.v, hr.v, col = colors, xlab = 'climate change', ylab = 'biodiversity risk', zlab = 'host richness'))

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")
pdf("CBDmaps.pdf",height = 10, width = 6)

clim.plot /
  bii.plot /
  disease.plot /
  risk.cm
dev.off()


