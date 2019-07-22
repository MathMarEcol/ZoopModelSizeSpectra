# Function which adds satellite SST, Chl and Bathy data to any data frame 
# with "Latitude" and "Longitude" variables
#
# Input:
# dat - dataframe with Lat/Lon variables.
# EnviroDir - directory where the MODIS SST and CHL files are
# bathy_file - link to GEBOC bathymetry nc file
#
# Jason D Everett (Mon 1st July 2019)

library(raster)
fAddEnviro <- function(dat, EnviroDir, bathy_file){
  
  dat <- dat %>% add_column(Bathy = NA, SST = NA, Chl = NA)
  
  ########### Import CHLOROPHYLL AND SST ###########
  SST_files <- list.files(path = EnviroDir, pattern = "MC_SST", recursive = FALSE, full.names = TRUE)
  Chl_files <- list.files(path = EnviroDir, pattern = "MC_CHL", recursive = FALSE, full.names = TRUE)
  
  # Ensure they are in order
  SST_files <- sort(SST_files)
  Chl_files <- sort(Chl_files)
  
  pb <- txtProgressBar(min = 0, max = 12, initial = 0, style = 3)
  
  for(i in 1:12){
    
    datmth <- dat %>% 
      filter(Month==i)
    
    rasSST <- raster(SST_files[i], varname = "sst", nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90)
    sst <- extract(rasSST, datmth[,c("Longitude", "Latitude")], na.rm=T, method = "simple")
    datmth$SST <- as.numeric(sst)
    
    rasChl <- raster(Chl_files[i], varname = "chlor_a", nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90)
    chl <- extract(rasChl, datmth[,c("Longitude", "Latitude")], na.rm=T, method = "simple")
    datmth$Chl <- as.numeric(chl)
    
    dat <- dat %>% 
      filter(Month!=i) %>% # Remove the old data from i month
      rbind(datmth) # Add in new data from i month
    
    setTxtProgressBar(pb, i)
    
    rm(rasSST, rasChl, sst, chl)
  }
  
  rasBathy <- raster(bathy_file, varname = "elevation")
  bathy <- extract(rasBathy, dat[,c("Longitude", "Latitude")], na.rm=T, method = "simple")
  dat$Bathy <- bathy
  
  return(dat)
  
}
