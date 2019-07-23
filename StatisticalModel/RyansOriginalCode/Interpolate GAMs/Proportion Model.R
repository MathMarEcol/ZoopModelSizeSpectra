rm(list=ls())
library(mgcv)
library(effects)
library(splines)
library(visreg)
library(raster)
library(maps)
library(mapdata)
library(maptools)
library(RColorBrewer)
library(stringr)
library(plyr)
library(RColorBrewer)
library(rgdal)
library("mgcv")

work_direct_main = "~/Library/Mobile Documents/com~apple~CloudDocs/PhD/" # Set your main directory where you
# put "Zooplankton GAMs"

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Zooplankton GAMs/Interpolate GAMs")
source("Zoo Prediction Function.R")


####################################################
############# IMPORT GAMS ##########################
chaetogs <- readRDS("chaets_gam11.rds") # Chaetognaths
carncops <- readRDS("copecarn_gam11.rds") # Carn Copepods
omnicops <- readRDS("omnicops_gam11.rds") # Omni/Herb Copepods 
totcops <- readRDS("copetotal_gam11.rds") # Total Copepods
krill <- readRDS("krill_gam11.rds") # Krill
salps <- readRDS("salps_gam11.rds") # Salps
larvs <- readRDS("larvs_gam11.rds") # Larvaceans
jelly <- readRDS("jelly_gam11.rds") # Jellyfish

####################################################
############# GLOBAL DATA ###########################
bathy <- read.csv("bathy.csv")

sst_clim <- read.csv("sst_clim_2.csv")
year_ave_sst = rowMeans(sst_clim[,c(3:14)], na.rm = TRUE)

chlo_clim <- read.csv("chlo_clim.csv")
year_ave_chlo = rowMeans(chlo_clim[,c(3:14)], na.rm = TRUE)

####### Create frames of zeros for bathymetry, sst and chlo - if you're wanting
####### to hold them constant. Note if you run sst constant, the model_zoo function
####### also holds day of year constant at 1

bathy_zero <- bathy
bathy_zero[,3] <- 0

sst_zero <- sst_clim
sst_zero[,c(3:14)] <- 0

chlo_zero <- chlo_clim
chlo_zero[,c(3:14)] <- 0

####################################################
###### RUN THE GAMs ################################

## You can run GEAR_MESH as either 'max' or 'mean'. Max finds the abundances using the gear_mesh level 
## with the largest abundances, and mean calculates the average abundance across all gear_mesh levels

krillt <- model_zoo(curr_gam = krill, BATH = bathy, sst_c = sst_clim,
                    chlo_c = chlo_clim, GEAR_MESH = 'mean') 

chaetot <- model_zoo(curr_gam = chaetogs, BATH = bathy, sst_c = sst_clim,
                     chlo_c = chlo_clim, GEAR_MESH = 'mean')

jellyt <- model_zoo(curr_gam = jelly, BATH = bathy, sst_c = sst_clim,
                    chlo_c = chlo_clim, GEAR_MESH = 'mean') 

salp_t <-  model_zoo(curr_gam = salps, BATH = bathy, sst_c = sst_clim,
                     chlo_c = chlo_clim, GEAR_MESH = 'mean') 

larv_t <-  model_zoo(curr_gam = larvs, BATH = bathy, sst_c = sst_clim,
                     chlo_c = chlo_clim, GEAR_MESH = "mean") 

omnic_t <- model_zoo(curr_gam = omnicops, BATH = bathy, sst_c = sst_clim,
                     chlo_c = chlo_clim, GEAR_MESH = "mean") 

carnc_t <- model_zoo(curr_gam = carncops, BATH = bathy, sst_c = sst_clim,
                     chlo_c = chlo_clim, GEAR_MESH = "mean") 

#cope_tot_t <- model_zoo(curr_gam = totcops, BATH = bathy, sst_c = sst_clim,
#                        chlo_c = chlo_clim, GEAR_MESH = 'max') 

## YEAR AVE ABUNDANCES
chaetss <- 10^chaetot$year_ave
 
krillss <- 10^krillt$year_ave

jellyss <- 10^jellyt$year_ave

larvss <- 10^larv_t$year_ave

salpss <- 10^salp_t$year_ave

ocopess <- 10^omnic_t$year_ave

ccopess <- 10^carnc_t$year_ave

#tcopess <- 10^cope_tot_t$year_ave

tot_zoo = chaetss + krillss + jellyss + larvss + salpss + ocopess + ccopess


## PLOT GLOBAL MAPS OF ABUNDANCES AND PROPORTIONS
lat_ranges = data.frame("min" = c(-38, -48, -43, -50, -42, -38, -40), "max" = c(38, 48, 43, 50, 42, 38, 40))

abund_data = data.frame("lat" = bathy$lat, "lon" = bathy$lon, "Larvaceans" = larvss, 
                        "Omni Cop" = ocopess, "Carn Cop" = ccopess,
                        "Chaetog" = chaetss, "Euphausiids" = krillss, "Salps" = salpss,
                        "Jellyfish" = jellyss)

plotss <- plot_all_abunds(abund_data, lat_ranges, "MEAN_GEAR_MAPS")
ggsave(filename = paste("ABUND_PROP_MAPS_MEAN_GEAR", "1.png", sep = ""), plot = ggarrange(plots = plotss[1:8], nrow = 4), width = 9, height = 9)
ggsave(filename = paste("ABUND_PROP_MAPS_MEAN_GEAR", "2.png", sep = ""), plot = ggarrange(plots = plotss[9:14], nrow = 3), width = 9, height = 7)

## PLOT PROPORTIONS AGAINST SST AND CHLO, USEFUL IF YOU ARE LOOKING AT MAIN EFFECTS
prop_frame = data.frame( "Larvs" = larvss/tot_zoo,"Omni_Cops" = ocopess/tot_zoo, "Carn_Cops" = ccopess/tot_zoo, 
                        "Chaets" = chaetss/tot_zoo,  "Krill" = krillss/tot_zoo, "Salps" = salpss/tot_zoo,
                        "Jellys" = jellyss/tot_zoo)

zoo_names <- c("Larvs","Omni_Cops", "Carn_Cops",  "Chaets", "Krill","Salps","Jellys")

plot_name = paste('Zoo_Stat_Proportions', "_chlo.png")
png(plot_name, width = 12, height = 16, units = "in", res = 72)
par(mfrow = c(4,2))
par(mfrow = c(4,2))
for(i in 1:length(zoo_names)){
  plot(log10(year_ave_chlo), prop_frame[,i], 
       xlab = 'log10(Chlorophyll)', ylab = 'Proportion', main = zoo_names[i])
}
dev.off()

plot_name = paste('Zoo_Stat_Proportions', "_sst.png")
png(plot_name, width = 12, height = 16, units = "in", res = 72)
par(mfrow = c(4,2))
for(i in 1:length(zoo_names)){
  plot(year_ave_sst, prop_frame[,i], 
       xlab = 'SST', ylab = 'Proportion', main = zoo_names[i])
}
dev.off()

par(mfrow = c(1,1))


############# TOTAL GLOBAL PIE CHART AND CENSUS #######################
slices <- c(sum(ocopess, na.rm = TRUE), 
            sum(ccopess, na.rm = TRUE),
           sum(jellyss, na.rm = TRUE),
            sum(larvss, na.rm = TRUE),
            sum(krillss, na.rm = TRUE),
           sum(chaetss, na.rm = TRUE),
            sum(salpss, na.rm = TRUE))
lbls <- c("O Cop", "C Cop", "Jellys", "Larvs", "Krill", "Chaets","Salps")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
   main = "Top 200m Global Zoo Abund Mix (> 200um)", radius = 1) 

############# OLIGO #######################
oligo <- which(log10(year_ave_chlo) <= 0.1)
slices <- c(sum(ocopess[oligo], na.rm = TRUE), 
            sum(ccopess[oligo], na.rm = TRUE),
            sum(jellyss[oligo], na.rm = TRUE),
            sum(larvss[oligo], na.rm = TRUE),
            sum(krillss[oligo], na.rm = TRUE),
            sum(chaetss[oligo], na.rm = TRUE),
            sum(salpss[oligo], na.rm = TRUE))
lbls <- c("O Cop", "C Cop", "Jellys", "Larvs", "Krill", "Chaets","Salps")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main = "Top 200m Oligo Global Zoo \n Abund Mix (> 200um and Chlo <= 0.1)", radius = 1) 


############# EUTRO REGIONS #######################
eutro <- which(log10(year_ave_chlo) > 1)
slices <- c(sum(ocopess[eutro], na.rm = TRUE), 
            sum(ccopess[eutro], na.rm = TRUE),
            sum(jellyss[eutro], na.rm = TRUE),
            sum(larvss[eutro], na.rm = TRUE),
            sum(krillss[eutro], na.rm = TRUE),
            sum(chaetss[eutro], na.rm = TRUE),
            sum(salpss[eutro], na.rm = TRUE))
lbls <- c("O Cop", "C Cop", "Jellys", "Larvs", "Krill", "Chaets","Salps")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main = "Top 200m Eutro Global Zoo \n Abund Mix (> 200um and Chlo >= 1)", radius = 1) 

#### TOT ZOO
library(raster)
r <- raster()
areas <- as.data.frame(area(r))*1e6
tot_ocops <- sum((ocopess*200)*areas/10, na.rm = TRUE)
tot_ccops <- sum((ccopess*200)*areas/10, na.rm = TRUE)
tot_salps <- sum((salpss*200)*areas/10, na.rm = TRUE)
tot_larvs <- sum((larvss*200)*areas/10, na.rm = TRUE)
tot_krill<- sum((krillss*200)*areas/10, na.rm = TRUE)
tot_jelly <- sum((jellyss*200)*areas/10, na.rm = TRUE)
tot_chaet <- sum((chaetss*200)*areas/10, na.rm = TRUE)

tot_zoo <- tot_ocops + tot_ccops + tot_salps + tot_larvs +
  tot_krill + tot_jelly + tot_chaet


###### AGGREGATE TO 5X5 GRID SQUARES, FOR COMPARISON WITH SIZE SPECTRUM MODEL
# 5x5 degree grid squares globally
lon <- seq(-177.5,177.5,5)
lat <- seq(-87.5,87.5,5)
lonlat <- expand.grid(lon = lon, lat = lat)

coordinates(lonlat) <- ~lon+lat
crs(lonlat) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

grd_lonlat <- points2grid(lonlat)
grd_lonlat <- SpatialGrid(grd_lonlat, proj4string = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

### CHAETOGNATHS
chaet5 <- data.frame("lon" = chaetot$lon, "lat" = chaetot$lat, "year_ave" = chaetot$year_ave)
coordinates(chaet5) <- ~ lon + lat
crs(chaet5) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
chaet_5 <- raster(aggregate(chaet5, grd_lonlat, mean, na.rm = TRUE))
chaet_5d <- as.data.frame(chaet_5, xy = TRUE)

### KRILL
krill5 <- data.frame("lon" = krillt$lon, "lat" = krillt$lat, "year_ave" = krillt$year_ave)
coordinates(krill5) <- ~ lon + lat
crs(krill5) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
krill_5 <- raster(aggregate(krill5, grd_lonlat, mean, na.rm = TRUE))
krill_5d <- as.data.frame(krill_5, xy = TRUE)

### JELLYFISH
jelly5 <- data.frame("lon" = jellyt$lon, "lat" = jellyt$lat, "year_ave" = jellyt$year_ave)
coordinates(jelly5) <- ~ lon + lat
crs(jelly5) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
jelly_5 <- raster(aggregate(jelly5, grd_lonlat, mean, na.rm = TRUE))
jelly_5d <- as.data.frame(jelly_5, xy = TRUE)

### SALPS
salp5 <- data.frame("lon" = salp_t$lon, "lat" = salp_t$lat, "year_ave" = salp_t$year_ave)
coordinates(salp5) <- ~ lon + lat
crs(salp5) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
salp5 <- raster(aggregate(salp5, grd_lonlat, mean, na.rm = TRUE))
salp_5d <- as.data.frame(salp5, xy = TRUE)

### LARVACEANS
larv5 <- data.frame("lon" = larv_t$lon, "lat" = larv_t$lat, "year_ave" = larv_t$year_ave)
coordinates(larv5) <- ~ lon + lat
crs(larv5) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
larv5 <- raster(aggregate(larv5, grd_lonlat, mean, na.rm = TRUE))
larv_5d <- as.data.frame(larv5, xy = TRUE)

### OMNI COPS
ocop5 <- data.frame("lon" = omnic_t$lon, "lat" = omnic_t$lat, "year_ave" = omnic_t$year_ave)
coordinates(ocop5) <- ~ lon + lat
crs(ocop5) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ocop5 <- raster(aggregate(ocop5, grd_lonlat, mean, na.rm = TRUE))
ocop_5d <- as.data.frame(ocop5, xy = TRUE)

### CARN COPS
ccop5 <- data.frame("lon" = carnc_t$lon, "lat" = carnc_t$lat, "year_ave" = carnc_t$year_ave)
coordinates(ccop5) <- ~ lon + lat
crs(ccop5) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ccop5 <- raster(aggregate(ccop5, grd_lonlat, mean, na.rm = TRUE))
ccop_5d <- as.data.frame(ccop5, xy = TRUE)

year_ave_data <- data.frame("Lon" = ccop_5d$x, "Lat" = ccop_5d$y, "Larvaceans" = larv_5d$year_ave, 
                            "Omni.Cop" = ocop_5d$year_ave, "Carn.Cop" = ccop_5d$year_ave, "Krill" = krill_5d$year_ave,
                            "Chaetog" = chaet_5d$year_ave, "Salps" = salp_5d$year_ave, "Jellys" = jelly_5d$year_ave)
write.csv(year_ave_data, "5deg_year_ave_data.csv", row.names = FALSE)
