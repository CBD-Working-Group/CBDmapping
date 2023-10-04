#LOAD DATA

library(reshape2)
library(raster)
library(ggplot2)
library(viridisLite)
library(patchwork)
library(colormap)
library(fields)
library(rgl)
library(rgdal)

library(maptools)
data(wrld_simpl)


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
#setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/climate/MiaoHotspot/")
library(ncdf4)
chot <- raster("MiaoHotspot/SED_ssp245_2080-2099.nc")
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
#bii.r= readRDS("bii.rds")

#https://figshare.com/articles/dataset/Global_maps_of_Biodiversity_Intactness_Index_Sanchez-Ortiz_et_al_2019_-_bioRxiv_/7951415/1
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/biodiversity/BII2019/abundance_richness_BII_maps/")
bii2<- raster("~/Desktop/final-rich-bii-isl-main.tif")

bii.r = crop(bii2, extent.r)
bii.r<- aggregate(bii.r, fact=4)
rm(bii2)

#------------------
#compare patterns

#extract all to locate data
chot.v= raster::extract(chot.ll, c(1:ncell(chot.ll)))
inds= which(!is.na(chot.v))
chot.v= chot.v[inds]

xys= xyFromCell(chot.ll, inds)

#extract values
bii.v= raster::extract(bii.r, xys)
hr.v= raster::extract(hr.ll, xys)

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

library(tidyverse)
#rescale between 0 and 1 for additive map
xy.dat %>%
  mutate(chot = (chot.v-min(chot.v))/(max(chot.v)-min(chot.v)),
         bii.threat = (bii.threat-min(bii.threat))/(max(bii.threat)-min(bii.threat)),
         hr = (hr.v-min(hr.v))/(max(hr.v)-min(hr.v)), 
         add = chot + bii.threat + hr,
         scale = (add-min(add))/(max(add)-min(add))
         ) -> xy.dat.norm1

#additive plot
#world outline
use<-fortify(worldcrop, region="ISO2")
addplot<-ggplot(xy.dat.norm1, aes(x= x, y=y, fill=scale, color=scale))+
  geom_tile()+
  #scale_fill_identity()+
  scale_color_identity()+
  scale_fill_gradient(low = "#FFFFFF", high = "#381393",space = "Lab", na.value = "white",guide = "colourbar",
                      aesthetics = "fill")+
  theme_void()+
  geom_polygon(data=use, mapping=aes(x=long,y=lat, group=group), fill="transparent", color="black")+
  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"), legend.position = "bottom")+ 
  scale_x_continuous(limits=c(-170,180), expand=c(0,0)) +
  scale_y_continuous(limits=c(-55,75), expand=c(0,0)) +
  labs(x = NULL, y = NULL)

#rescale between 1 and 2 for multiplicative map
xy.dat %>%
  mutate(chot = (chot.v / max(chot.v)+1),
         bii.threat = (bii.threat / max(bii.threat)+1),
         hr = (hr.v / max(hr.v)+1), 
         mult = chot * bii.threat * hr,
         scale = (mult-min(mult))/(max(mult)-min(mult))
  ) -> xy.dat.norm2


#multiplicative plot
multplot <- ggplot(xy.dat.norm2, aes(x= x, y=y, color=scale, fill = scale))+
  geom_tile()+
  #scale_fill_identity()+
  scale_color_identity()+
  scale_fill_gradient(low = "#FFFFFF", high = "#381393",space = "Lab", na.value = "white",guide = "colourbar",
                      aesthetics = "fill")+
  theme_void()+
  geom_polygon(data=use, mapping=aes(x=long,y=lat, group=group), fill="transparent", color="black")+
  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"), legend.position = "bottom")+ 
  scale_x_continuous(limits=c(-170,180), expand=c(0,0)) +
  scale_y_continuous(limits=c(-55,75), expand=c(0,0)) +
  labs(x = NULL, y = NULL)


#look at the top 20% percent 
#quantile(xy.dat.norm1$chot, probs = .8) #0.6436111
#quantile(xy.dat.norm1$bii.threat, probs = .8) #0.6508289  
#quantile(xy.dat.norm1$hr, probs = .8) #0.4861111

#make new dataframe
new<-data.frame(x= xy.dat$x,
                y= xy.dat$y, 
                chot = NA, 
                bii = NA, 
                hr = NA)

for (i in 1:nrow(new)){
  if (xy.dat.norm1$chot[i] >= 0.6436111) {new$chot[i] = 1}
  if (xy.dat.norm1$bii.threat[i] >= 0.6508289) {new$bii[i] = 1}
  if (xy.dat.norm1$hr[i] >= 0.4861111) {new$hr[i] = 1}
}
#anything not in the top 20% is assigned a zero
new[is.na(new)] <- 0

xy.dat.perc<-new
xy.dat.perc$colors<-colors3d(xy.dat.perc[,c(3,4,5)])

#change colors to match the paper
xy.dat.perc$colors<-gsub("#FF0000", "#8C212A", xy.dat.perc$colors) #C
xy.dat.perc$colors<-gsub("#000000", "white",  xy.dat.perc$colors) #NONE
xy.dat.perc$colors<-gsub("#00FF00", "#214B8D", xy.dat.perc$colors) #B
xy.dat.perc$colors<-gsub("#0000FF", "#C7AB80", xy.dat.perc$colors) #ID
xy.dat.perc$colors<-gsub("#00FFFF", "#AAC5E2", xy.dat.perc$colors) #B + ID
xy.dat.perc$colors<-gsub("#FFFF00", "#53998A", xy.dat.perc$colors) #C + B
xy.dat.perc$colors<-gsub("#FF00FF", "#D78244", xy.dat.perc$colors) #C + ID
xy.dat.perc$colors<-gsub("#FFFFFF", "#381393", xy.dat.perc$colors) #CBD

#percentile plot
threshold <- ggplot(xy.dat.perc, aes(x=x, y=y, color=colors, fill =colors))+
  geom_tile()+
  scale_fill_identity()+
  scale_color_identity()+
  theme_void()+
  geom_polygon(data=use, mapping=aes(x=long,y=lat, group=group), fill="transparent", color="black")+
  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"), legend.position = "left")+ 
  scale_x_continuous(limits=c(-170,180), expand=c(0,0)) +
  scale_y_continuous(limits=c(-55,75), expand=c(0,0)) +
  labs(x = NULL, y = NULL)


layout <- "
CC
CC
AB
"
  
addplot + multplot +threshold + plot_layout(design = layout, guides = 'collect')
ggsave("map_figure.png", height=8, width =12)

#-----
#plot together
use<-fortify(worldcrop, region="ISO2")
#climate
chot.df <- as.data.frame(chot, xy=TRUE) #Convert raster to data.frame
names(chot.df)[3] <- 'climate' 
head(chot.df)

clim.plot= ggplot(data = chot.df)+
  geom_raster(mapping=aes(x=x, y=y, fill=climate))+
  #scale_fill_gradientn(colours= rev(terrain.colors(10)), name='climate change', na.value="white")+
  scale_fill_gradient(low = "#FFFFFF", high = "#990F26",space = "Lab", na.value = "white",guide = "colourbar",
    aesthetics = "fill")+
  geom_polygon(data=use, mapping=aes(x=long,y=lat, group=group), fill="transparent", color="black")+
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
  #scale_fill_gradientn(colours= rev(terrain.colors(10)), name='biodiversity risk', na.value="white")+
  scale_fill_gradient(low = "#FFFFFF", high = "#084C92",space = "Lab", na.value = "white",guide = "colourbar",
                      aesthetics = "fill")+
  geom_polygon(data=use, mapping=aes(x=long,y=lat, group=group), fill="transparent", color="black")+
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
  #scale_fill_gradientn(colours= rev(terrain.colors(10)), name='host richness', na.value="white")+ #'predicted host-parasite interactions'
  scale_fill_gradient(low = "#FFFFFF", high = "#99600F",space = "Lab", na.value = "white",guide = "colourbar",
                      aesthetics = "fill")+
  geom_polygon(data=use, mapping=aes(x=long,y=lat, group=group), fill="transparent", color="black")+
  theme_void()+
  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"), legend.position = "bottom")+ 
  scale_x_continuous(limits=c(-170,180), expand=c(0,0)) +
  scale_y_continuous(limits=c(-55,75), expand=c(0,0)) +
  labs(x = NULL, y = NULL)

#Risk map
# map color to the climate variables
xy.dat$colors <- colors3d(xy.dat[,c(3,4,5)]) #colors is additive
#names(xy.dat)[3:5]=c("climate change","biodiversity risk","host richness")

#plot
risk.cm=ggplot(xy.dat, aes(x= x, y=y, fill=colors, color=colors))+
  geom_tile()+
  scale_fill_identity()+scale_color_identity()+theme_void()+
  geom_polygon(data=use, mapping=aes(x=long,y=lat, group=group), fill="transparent", color="black")+
  theme(plot.margin=grid::unit(c(0,0,0,0), "cm"), legend.position = "bottom")+ 
  scale_x_continuous(limits=c(-170,180), expand=c(0,0)) +
  scale_y_continuous(limits=c(-55,75), expand=c(0,0)) +
  labs(x = NULL, y = NULL)

#-----
#plot 3D legend
with(xy.dat, plot3d(chot.v, bii.v, hr.v, theta = 20, phi = 20, col = colors, xlab = 'climate change', ylab = 'biodiversity risk', zlab = 'host richness'))

rgl.snapshot('3dplot.png', fmt = 'png')
#rgl.postscript('3dplot.pdf', fmt = 'pdf')



library(magick)
library(ggpubr)
leg3d <- image_read('3dplot.png')
leg3d <- ggplot() +
  background_image(leg3d) + coord_fixed()

#-----
#layout <- "
#AADDDD
#BBDDDD
#CC####
#"

layout <- "
AA####
AADDDD
BBDDDD
BBDDDD
CCDDDD
CCEE##
"
maps= clim.plot + bii.plot + disease.plot + risk.cm + leg3d + plot_layout(design = layout)
ggsave("figure2.png", height=8, width =12)

layout2 <- "
AABBCC
DDDDDD
DDDDDD
DDDDDD
"
maps2= clim.plot + bii.plot + disease.plot + risk.cm + plot_layout(design = layout2)
ggsave("figure2_l2.png", height=8, width =12)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/figs/")
pdf("CBDmaps.pdf",height = 6, width = 12)

maps + inset_element(leg3d, left = 1, bottom = 0, right = 0.6, top = 0.5, align_to = 'full')


dev.off()
