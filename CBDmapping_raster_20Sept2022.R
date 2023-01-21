#LOAD DATA

library(reshape2)
library(raster)
library(ggplot2)
library(viridisLite)
library(patchwork)
library(colormap)
library(fields)
library(rgl)

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
hr <- raster("Total host richness.grd")
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

xy.dat$bii.threat= 1-xy.dat$bii.v

#=========================
#3way color map
#http://matthewkling.github.io/media/colormap/
#devtools::install_github("matthewkling/colormap")

xy.dat= as.data.frame(cbind(xys, chot.v, bii.v, hr.v)) 
xy.dat$bii.threat= 1-xy.dat$bii.v
xy.dat= na.omit(xy.dat)

#write out
setwd(hdir)
#write.csv(xy.dat,"threatmapdata.csv")
#xy.dat= read.csv("threatmapdata.csv")

#-----
#plot together

#climate
chot.df <- as.data.frame(chot, xy=TRUE) #Convert raster to data.frame
names(chot.df)[3] <- 'climate' 
head(chot.df)

clim.plot= ggplot(data = chot.df)+
  geom_raster(mapping=aes(x=x, y=y, fill=climate))+
  scale_fill_gradientn(colours= rev(terrain.colors(10)), name='climate change', na.value="white")+
  theme_void()+
  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"), legend.position = "bottom")+ 
  scale_x_continuous(limits=c(-170,180), expand=c(0,0)) +
  scale_y_continuous(limits=c(-55,75), expand=c(0,0)) +
  labs(x = NULL, y = NULL)

#biodiversity
bii.risk= 1-bii.r
bii.df <- as.data.frame(bii.risk, xy=TRUE) #Convert raster to data.frame
names(bii.df)[3] <- 'biodiversity' 

bii.plot= ggplot(data = bii.df)+
  geom_raster(mapping=aes(x=x, y=y, fill=biodiversity))+
  scale_fill_gradientn(colours= rev(terrain.colors(10)), name='biodiversity risk', na.value="white")+
  theme_void()+
  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"), legend.position = "bottom")+ 
  scale_x_continuous(limits=c(-170,180), expand=c(0,0)) +
  scale_y_continuous(limits=c(-55,75), expand=c(0,0)) +
  labs(x = NULL, y = NULL)

#disease
hpar.df <- as.data.frame(hr.ll, xy=TRUE) #Convert raster to data.frame
names(hpar.df)[3] <- 'disease' 
head(hpar.df)

disease.plot= ggplot(data = hpar.df)+
  geom_raster(mapping=aes(x=x, y=y, fill=disease))+
  scale_fill_gradientn(colours= rev(terrain.colors(10)), name='host richness', na.value="white")+ #'predicted host-parasite interactions'
  theme_void()+
  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"), legend.position = "bottom")+ 
  scale_x_continuous(limits=c(-170,180), expand=c(0,0)) +
  scale_y_continuous(limits=c(-55,75), expand=c(0,0)) +
  labs(x = NULL, y = NULL)

#Risk map
# map color to the climate variables
xy.dat$colors <- colors3d(xy.dat[,c(3,4,5)])
#names(xy.dat)[3:5]=c("climate change","biodiversity risk","host richness")

risk.cm=ggplot(xy.dat, aes(x= x, y=y, fill=colors, color=colors))+geom_tile()+
  scale_fill_identity()+scale_color_identity()+theme_void()+
  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"), legend.position = "bottom")+ 
  scale_x_continuous(limits=c(-170,180), expand=c(0,0)) +
  scale_y_continuous(limits=c(-55,75), expand=c(0,0)) +
  labs(x = NULL, y = NULL)

#-----
#plot 3D legend
with(xy.dat, plot3d(chot.v, bii.v, hr.v, col = colors, xlab = 'climate change', ylab = 'biodiversity risk', zlab = 'host richness'))

rgl.snapshot('3dplot.png', fmt = 'png')
#rgl.postscript('3dplot.pdf', fmt = 'pdf')

library(magick)
library(ggpubr)
leg3d <- image_read('3dplot.png')
leg3d <- ggplot() +
  background_image(leg3d) + coord_fixed()

#-----
layout <- "
AADDD
BBDDD
CCDDD
"

maps= clim.plot + bii.plot + disease.plot + risk.cm +
  plot_layout(design = layout)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")
pdf("CBDmaps.pdf",height = 6, width = 12)

maps + inset_element(leg3d, left = 0.6, bottom = 0, right = 1, top = 0.7, align_to = 'full')

dev.off()
