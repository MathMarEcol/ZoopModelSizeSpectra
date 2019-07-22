# library(sp)
# library(rgdal)
library(tidyverse)
library(ncdf4)
library(raster)

EnviroDir <- "Data/EnvironmentalData"

## AGGREGATE TO 1X1 GRID SQUARES

mnth <- c("01January", "02February", "03March", "04April", "05May", 
          "06June", "07July", "08August", "09September", "10October", "11November", "12December")

##### ##### Create the new resolution to match the SSM ##### #####
# 1 x 1 degree grid squares globally

lon <- seq(-179.5,179.5,1)
lat <- seq(-89.5,89.5,1)
Coord_ref <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
lonlat <- expand.grid(Long = lon, Lat = lat)

coordinates(lonlat) <- ~Long+Lat
proj4string(lonlat) <- CRS(Coord_ref)
gridded(lonlat) = TRUE
grd_lonlat <- raster(lonlat)


##### ##### Do Bathy first ##### #####
bathy_file <- paste0(EnviroDir,"/GEBCO_2014_2D.nc")

rasBathy <- raster(bathy_file, varname = "elevation")
Bathyagg <- resample(rasBathy, grd_lonlat, method='bilinear')
Bathyagg <- reclassify(Bathyagg, c(0,Inf,NA))
Bathyagg <- abs(Bathyagg)


names(Bathyagg) <- "Bathy"
saveRDS(Bathyagg, file = paste0(EnviroDir,"/Bathy_raster_1Deg.rds"))


##### ##### No do SST and Chl ##### #####
SST_files <- list.files(path = EnviroDir, pattern = "MC_SST", recursive = FALSE, full.names = TRUE)
Chl_files <- list.files(path = EnviroDir, pattern = "MC_CHL", recursive = FALSE, full.names = TRUE)

for(i in 1:12){
  
  ########### Import CHLOROPHYLL AND SST ###########
  rasSST <- raster(SST_files[i], varname = "sst", nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90)
  SSTagg <- resample(rasSST, grd_lonlat, method='bilinear')
  names(SSTagg) <- "SST"
  saveRDS(SSTagg, file = paste0(EnviroDir,"/SST_raster_",mnth[i],"_1Deg.rds"))
  
  if (i == 1){
    SSTbrick <- brick(SSTagg)} else{
      SSTbrick <- addLayer(SSTbrick,SSTagg)
    }

  rasChl <- raster(Chl_files[i], varname = "chlor_a", nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90)
  Chlagg <- resample(rasChl, grd_lonlat, method='bilinear')
  names(Chlagg) <- "Chl"
  saveRDS(Chlagg, file = paste0(EnviroDir,"/Chl_raster_",mnth[i],"_1Deg.rds"))
  
  if (i == 1){
    Chlbrick <- brick(Chlagg)} else{
      Chlbrick <- addLayer(Chlbrick,Chlagg)
    }
}


SSTbrick_mn  <- calc(SSTbrick, fun = mean)
names(SSTbrick_mn) <- "SST"
saveRDS(SSTbrick_mn, file = paste0(EnviroDir,"/SST_raster_00Mean_1Deg.rds"))


Chlbrick_mn  <- calc(Chlbrick, fun = mean)
names(Chlbrick_mn) <- "Chl"
saveRDS(Chlbrick_mn, file = paste0(EnviroDir,"/Chl_raster_00Mean_1Deg.rds"))

## I still need to create a raster brick for all satellite data, but not sure if it is needed so I will leave it for the moment
# colnames(SST_store) <- colnames(Chl_store) <- c("Long","Lat","January", "February", "March", "April", "May", 
                                                # "June", "July", "August", "September", "October", "November", "December")
