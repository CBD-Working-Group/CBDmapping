# HP prediction risk maps

# 1. Subset posterior of full data / full model HP to links that are undocumented
# 2. Weight interactions by probability then sum per host
# 3. Make raster and plot by host (like weighted richness map)

require(letsR)
require(RColorBrewer)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/CBDwg/data/disease/HostParasite/")

# load observed network
com <- readRDS("./clean_data/full_com.rds")

# Load Prediction Data
load("./model_results/FULL-FULLSET.RData")

# order P to match com
P <- P[rownames(com), colnames(com)]

# Creating posterior probability matrix with only undocumented interactions 
P_undoc <- P
P_undoc[com==1] <- 0

#----
# Read shape file
# Not included in data & code supplement because of size - download from IUCN RedList directly
shp <- rgdal::readOGR("./raw_data/IUCN/Mammals_Terrestrial.shp")
str(shp@data)

shp@data$BINOMIAL <- gsub(" ","_", as.character(shp@data$BINOMIAL))

mamms.in.DB <- which(unique(as.character(shp@data$BINOMIAL)) %in% rownames(P_undoc))
length(mamms.in.DB)
# 1674 species in dat with shape files

mamms.missing <- sort(unique(rownames(P_undoc)[which(!rownames(P_undoc) %in% as.character(shp@data$BINOMIAL))]))
length(mamms.missing)
# 161 missing host species

mamms.in.map <- sort(unique(as.character(shp@data$BINOMIAL)))
#----

# Trying to merge missing using the Phylacine v1.2 synonym table
syn <- read.csv("./raw_data/Synonymy_table_with_unaccepted_species.csv", as.is=T)
syn$Elton_binomial <- paste0(syn$EltonTraits.1.0.Genus,"_",syn$EltonTraits.1.0.Species)
intersect(mamms.missing,syn$Elton_binomial)# all 161

syn$IUCN_binomial <- paste0(syn$IUCN.2016.3.Genus,"_",syn$IUCN.2016.3.Species)
missing_iucn_names <- syn$IUCN_binomial[syn$Elton_binomial%in%mamms.missing]
intersect(missing_iucn_names, mamms.missing)#92

# Swapping names in com (easier as the map is just for illustration purposes)
lookup <- setNames(syn$IUCN_binomial,syn$Elton_binomial)
lookup[mamms.missing]

rownames(P_undoc)[rownames(P_undoc)%in%mamms.missing] <- lookup[rownames(P_undoc)[rownames(P_undoc)%in%mamms.missing]]

# Verifying it worked
mamms.in.DB <- which(unique(as.character(shp@data$BINOMIAL)) %in% rownames(P_undoc))
length(mamms.in.DB)
# now 1724 species in dat with shape files

mamms.missing <- sort(unique(rownames(P_undoc)[which(!rownames(P_undoc) %in% as.character(shp@data$BINOMIAL))]))
length(mamms.missing)
# now 97 missing host species

# removing fully aquatic mammals 
mamms.missing <- mamms.missing[grep(
  "Balaena|Balaenoptera|Berardius|Cephalorhynchus|Delphinapterus|Delphinus|Dugong|Eubalaena|Eschrichtius|Feresa|Globicephala|Grampus|Hyperoodon|Inia|Kogia|Lagenodelphis|Lagenorhynchus|Lissodelphis|Megaptera|Mesoplodon|Monodon|Neophocaena|Orcaella|Orcinus|Peponocephala|Phocoena|Phocoenoides|Physeter|Platanista|Pontoporia|Pseudorca|Sotalia|Sousa|Stenella|Steno|Trichechus|Tursiops|Ziphius", 
  mamms.missing, invert=TRUE, ignore.case=FALSE)]
head(mamms.missing,30)

# Removing species with obviously no terrestrial IUCN range maps
mamms.missing <- mamms.missing[grep(
  "accepted|Bos_primigenius|Camelus_dromedarius|Homo_sapiens", 
  mamms.missing, invert=TRUE, ignore.case=FALSE)]
head(mamms.missing,30)

# Swapping for IUCN accepted synonyms (evaluated using IUCN online portal July 3rd 2019)
rownames(P_undoc)[rownames(P_undoc)=="Catopuma_temminckii"] <- "Pardofelis_temminckii"
rownames(P_undoc)[rownames(P_undoc)=="Cercopithecus_denti"] <- "Cercopithecus_pogonias"
rownames(P_undoc)[rownames(P_undoc)=="Chaerephon_plicatus"] <- "Tadarida_plicata"
rownames(P_undoc)[rownames(P_undoc)=="Chaerephon_pumilus"] <- "Tadarida_pumila"
rownames(P_undoc)[rownames(P_undoc)=="Galagoides_demidoff"] <- "Galagoides_demidovii"
rownames(P_undoc)[rownames(P_undoc)=="Hydrictis_maculicollis"] <- "Lutra_maculicollis"
rownames(P_undoc)[rownames(P_undoc)=="Hypsugo_savii"] <- "Pipistrellus_savii"
rownames(P_undoc)[rownames(P_undoc)=="Lycalopex_culpaeus"] <- "Pseudalopex_culpaeus"
rownames(P_undoc)[rownames(P_undoc)=="Lycalopex_fulvipes"] <- "Pseudalopex_fulvipes"
rownames(P_undoc)[rownames(P_undoc)=="Lycalopex_griseus"] <- "Pseudalopex_griseus"
rownames(P_undoc)[rownames(P_undoc)=="Lycalopex_gymnocercus"] <- "Pseudalopex_gymnocercus"
rownames(P_undoc)[rownames(P_undoc)=="Lycalopex_vetulus"] <- "Pseudalopex_vetulus"
rownames(P_undoc)[rownames(P_undoc)=="Mico_chrysoleucos"] <- "Mico_chrysoleucus"
rownames(P_undoc)[rownames(P_undoc)=="Mops_condylurus"] <- "Tadarida_condylura"
rownames(P_undoc)[rownames(P_undoc)=="Neoromicia_helios"] <- "Pipistrellus_helios"
# Neovision macrodon is extinct...
rownames(P_undoc)[rownames(P_undoc)=="Piliocolobus_badius"] <- "Procolobus_badius"
rownames(P_undoc)[rownames(P_undoc)=="Piliocolobus_kirkii"] <- "Procolobus_kirkii"
rownames(P_undoc)[rownames(P_undoc)=="Piliocolobus_rufomitratus"] <- "Procolobus_rufomitratus"
# Piliocolobus_tephrosceles is debated to be a subspecies of a number of different colobus species
rownames(P_undoc)[rownames(P_undoc)=="Rhinolophus_hildebrandtii"] <- "Rhinolophus_hildebrandti"
rownames(P_undoc)[rownames(P_undoc)=="Vulpes_lagopus"] <- "Alopex_lagopus"

# Presence-Absence matrix
mamms.in.DB <- intersect(unique(as.character(shp@data$BINOMIAL)), rownames(P_undoc))
length(mamms.in.DB)
# 1743 species in with shape files

shp <- subset(shp, BINOMIAL%in%mamms.in.DB)



# Create host presence-absence map
if (!file.exists("../clean_data/PresAbs.rds")) {
  
  ## generating pres-abs matrix
  PresAbs <- lets.presab(shp, 
                         xmn = -180, xmx = 180,  
                         resol=0.5, count=T) 
  
  saveRDS(PresAbs, "../clean_data/PresAbs.rds")
  
} else { PresAbs <- readRDS("../clean_data/PresAbs.rds") }


##################################
## MAP OF OBSERVED INTERACTIONS ##
##################################

cols <- "PuBuGn"
pal <- c(rep("#FFFFFF",3), colorRampPalette(brewer.pal(8,cols))(67))

para_richness <- rowSums(com)
para_richness <- para_richness[names(para_richness)%in%PresAbs$Species_name] # remove missing species
para_richness <- para_richness[PresAbs$Species_name] # re-order to match PresAbs

PresAbs_com <- lets.subsetPAM(PresAbs, names(para_richness), remove.cells = FALSE)

obs_map <- lets.maplizer(PresAbs_com, para_richness, names(para_richness),
                         func = sum, 
                         ras = TRUE)

# scaled_obs_map <- obs_map$Raster/cellStats(obs_map$Raster, "max")

pdf("./plots_tables/observed_map.pdf", width=11, height=6)
plot(obs_map$Raster, 
     # col = colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(256),
     col = pal,
     # col = colorRampPalette(colors)(100),
     xlim = c(-180, 180), ylim = c(-66, 84),
     axes = FALSE,
     box=FALSE,
     main = "Map of documented host-parasite interactions")
dev.off()


####################
## COMBINED MODEL ##
####################
load("../model_results/FULL-FULLSET.RData")
# Creating posterior probability matrix with only undocumented interactions 
P <- P[rownames(com), colnames(com)]# order P to match com
P_undoc_comb <- P
P_undoc_comb[com==1] <- 0
rownames(P_undoc_comb) <- rownames(P_undoc)

# Getting risk of undocumented parasite interactions per host species
risk_comb <- rowSums(P_undoc_comb)
risk_comb <- risk_comb[names(risk_comb)%in%PresAbs$Species_name] # remove missing species
risk_comb <- risk_comb[PresAbs$Species_name] # re-order to match PresAbs

riskmap_comb <- lets.maplizer(PresAbs, risk_comb, names(risk_comb),
                              func = sum,
                              ras = TRUE)

# rescaling raster to 0-1
scaled_riskmap_comb <- riskmap_comb$Raster/cellStats(riskmap_comb$Raster, "max")

# pdf("riskmap_comblogenymodel.pdf", width=11, height=6)
plot(scaled_riskmap_comb, 
     # col = colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(256),
     col = colorRampPalette(brewer.pal(9,"BuPu"))(256),
     # col = colorRampPalette(colors)(100),
     xlim = c(-180, 180), ylim = c(-66, 84),
     axes = FALSE,
     box=FALSE,
     main = "Hotspots of undocumented host-parasite interations (Combined model)")
# dev.off()


#####################
## PHYLOGENY MODEL ##
#####################

load("../model_results/DISTANCE-FULLSET.RData")
# Creating posterior probability matrix with only undocumented interactions 
P <- P[rownames(com), colnames(com)]# order P to match com
P_undoc_phy <- P
P_undoc_phy[com==1] <- 0
rownames(P_undoc_phy) <- rownames(P_undoc)

# Getting risk of undocumented parasite interactions per host species
risk_phy <- rowSums(P_undoc_phy)
risk_phy <- risk_phy[names(risk_phy)%in%PresAbs$Species_name] # remove missing species
risk_phy <- risk_phy[PresAbs$Species_name] # re-order to match PresAbs

riskmap_phy <- lets.maplizer(PresAbs, risk_phy, names(risk_phy),
                             func = sum,
                             ras = TRUE)

# rescaling raster to 0-1
scaled_riskmap_phy <- riskmap_phy$Raster/cellStats(riskmap_phy$Raster, "max")

# pdf("riskmap_phylogenymodel.pdf", width=11, height=6)
plot(scaled_riskmap_phy, 
     # col = colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(256),
     col = colorRampPalette(brewer.pal(9,"BuPu"))(256),
     # col = colorRampPalette(colors)(100),
     xlim = c(-180, 180), ylim = c(-66, 84),
     axes = FALSE,
     box=FALSE,
     main = "Hotspots of undocumented host-parasite interations (Phylogeny model)")
# dev.off()




####################
## AFFINITY MODEL ##
####################

load("../model_results/AFFINITY-FULLSET.RData")
# Creating posterior probability matrix with only undocumented interactions 
P <- P[rownames(com), colnames(com)]# order P to match com
P_undoc_aff <- P
P_undoc_aff[com==1] <- 0
rownames(P_undoc_aff) <- rownames(P_undoc)

# Getting risk of undocumented parasite interactions per host species
risk_aff <- rowSums(P_undoc_aff)
risk_aff <- risk_aff[names(risk_aff)%in%PresAbs$Species_name] # remove missing species
risk_aff <- risk_aff[PresAbs$Species_name] # re-order to match PresAbs

riskmap_aff <- lets.maplizer(PresAbs, risk_aff, names(risk_aff),
                             func = sum,
                             ras = TRUE)

# rescaling raster to 0-1
scaled_riskmap_aff <- riskmap_aff$Raster/cellStats(riskmap_aff$Raster, "max")

# pdf("riskmap_afflogenymodel.pdf", width=11, height=6)
plot(scaled_riskmap_aff, 
     # col = colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(256),
     col = colorRampPalette(brewer.pal(9,"BuPu"))(256),
     # col = colorRampPalette(colors)(100),
     xlim = c(-180, 180), ylim = c(-66, 84),
     axes = FALSE,
     box=FALSE,
     main = "Hotspots of undocumented host-parasite interations (Affinity model)")
# dev.off()


###################
# COMBINED FIGURE #
###################

cols <- "PuBuGn"
pal <- c(rep("#FFFFFF",3), colorRampPalette(brewer.pal(8,cols))(67))

pdf("../plots_tables/riskmap_tryptych.pdf", width=5, height=5.5)

par(mfrow=c(3,1), mai = c(0, 0.1, 0, 0.1) )

plot(scaled_riskmap_aff, 
     col = pal,
     xlim = c(-180, 180), ylim = c(-55, 84),
     axes = FALSE,
     box=FALSE)

title("A) Affinity model", line = -1.2, cex.main=1.15, family="sans", adj=0)

plot(scaled_riskmap_phy, 
     col = pal,
     xlim = c(-180, 180), ylim = c(-55, 84),
     axes = FALSE,
     box=FALSE)

title("B) Phylogeny model", line = -1.2, cex.main=1.15, family="sans", adj=0)

plot(scaled_riskmap_comb, 
     col = pal,
     xlim = c(-180, 180), ylim = c(-55, 84),
     axes = FALSE,
     box=FALSE)

title("C) Combined model", line = -1.2, cex.main=1.15, family="sans", adj=0)

dev.off()