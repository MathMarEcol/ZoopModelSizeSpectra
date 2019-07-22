rm(list = ls())

library(tidyverse)
library(ncdf4)
library(yaImpute)
library(lubridate)

EnviroDir <- "Data/EnvironmentalData"

########### IMPORT BATHYMETRY ###########
bathy_files <- "GEBCO_BATHY_2002-01-01_rgb_360x180.csv"
bathy_raw <- read.csv(paste0(EnviroDir,"/",bathy_files), header = FALSE)
bathy_vector <- as.vector(t(as.matrix(bathy_raw))) # Convert to vector

# create a grid of lat/lon, add headers, and convert to data frame
bathy_matrix = matrix(NA, nrow = length(bathy_vector), ncol = 3)
bathy_matrix[,c(1,2)] = as.matrix(expand.grid("Long" = seq(-179.5, 179.5, 1), "Lat" = seq(89.5, -89.5, -1)))
bathy_matrix[,3] <- bathy_vector
colnames(bathy_matrix) <- c("Long", "Lat", "BATHY")
bathy_matrix <- as.data.frame(bathy_matrix)
bathy_matrix[c(bathy_matrix$BATHY > 0), "BATHY"] <- NA
bathy_matrix <- bathy_matrix[!is.na(bathy_matrix$BATHY),] # Remove all the land so the nearest neighbbour below always finds a value
########### END BATHYMETRY ###########

########### Import CHLOROPHYLL AND SST ###########
SST_files <- list.files(path = EnviroDir, pattern = "MC_SST", recursive = TRUE, full.names = TRUE)
Chl_files <- list.files(path = EnviroDir, pattern = "MC_CHL", recursive = TRUE, full.names = TRUE)

# Rearrange to get January first up. 
SST_files <- SST_files[c(7:12, 1:6)]
Chl_files <- Chl_files[c(7:12, 1:6)]

# IMPORT NETCDF AND EXTRACT GEOGRAPHIC REFERENCES
nc.SST <- nc_open(SST_files[1])

SST_lat <- ncvar_get(nc.SST, "lat")
SST_lon <- ncvar_get(nc.SST, "lon")
SST_var <- t(ncvar_get(nc.SST, "sst4"))

# Create a matrix of data
SST_store <- Chl_store <- matrix(NA, nrow = length(SST_lat)*length(SST_lon), ncol = 14)
SST_store[,c(1,2)] <- Chl_store[,c(1,2)] <- as.matrix(expand.grid(lon = SST_lon, lat = SST_lat))

### OPEN NETCDFS
SST_nc_all <- lapply(SST_files, nc_open)
Chl_nc_all <- lapply(Chl_files, nc_open)

## EXTRACT SST AND CHLOROPHYLL, THESE SST AND CHLO CLIMATOLOGIES ARE 9KM MODIS AQUA
SST_nc_vall <- lapply(SST_nc_all, ncvar_get, varid = "sst4")
Chl_nc_vall <- lapply(Chl_nc_all, ncvar_get, varid = "chlor_a")


## PUT MONTHLY SST AND CHLO IN RESPECTIVE MATRICES, LABEL COLUMNS AND CONVERT TO DATA FRAME,
for(i in 1:12){
  SST_store[,c(i+2)] <- as.vector((unlist(SST_nc_vall[[i]])))
  Chl_store[,c(i+2)] <- as.vector((unlist(Chl_nc_vall[[i]])))
}
colnames(SST_store) <- colnames(Chl_store) <- c("Long","Lat","January", "February", "March", "April", "May", 
                                                "June", "July", "August", "September", "October", "November", "December")
SST_store <- as.data.frame(SST_store)
Chl_store <- as.data.frame(Chl_store)
########### End CHLOROPHYLL AND SST ###########


########### IMPORT ZOOPLANKTON DATA ###########
dat <- readRDS("LatestDatabaseOuput_Final.rds")
dat$Day[is.na(dat$Day)] <- 15 # Make all missing days to the middle of the month
dat <- dat[!is.na(dat$Month),] # Remove rows with no complete date

dat$Month2 <- month(dat$Month,label = TRUE, abbr = FALSE)

dat$SST = NA
dat$Chl = NA
dat$Bathy = NA
nvsc <- ann(as.matrix(SST_store[,c("Long", "Lat")]), as.matrix(dat[,c("Longitude", "Latitude")]), k = 1, verbose = FALSE)$knnIndexDist[,1]
nvb <- ann(as.matrix(bathy_matrix[,c("Long", "Lat")]), as.matrix(dat[,c("Longitude", "Latitude")]), k = 1, verbose = FALSE)$knnIndexDist[,1]

pb <- txtProgressBar(min = 0, max = dim(dat)[1], style = 3)
for(j in 1:dim(dat)[1]){
  setTxtProgressBar(pb, j)
  dat$SST[j] <- SST_store[nvsc[j], dat$Month2[j]]
  dat$Chl[j] <- Chl_store[nvsc[j], dat$Month2[j]]
  dat$Bathy[j] <- bathy_matrix[nvb[j],"BATHY"]
}

saveRDS(dat, file = "LatestDatabaseOuput_Final_Enviro.rds")
