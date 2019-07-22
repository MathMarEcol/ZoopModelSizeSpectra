### Author: Ryan Heneghan
### Updated by Jason Everett
### Date: May 2019
### Last Updated: 4th July 2019
### This script creates figures comparing output from the size spectrum model, with
### output from the statistical models for different zooplankton groups.

library(raster)
library(tidyverse)

############# ############# ############# ############# ############# 
#############   PREPARE DATA FROM NUMERICAL MODEL       #############
############# ############# ############# ############# ############# 
load("../../NumericalModel/Multi-ZooSizeSpectrumModel/31MayModelOutputFull/may30_raw_full.RData")
ssm <- may30_raw_full
rm(may30_raw_full)
# 'data.frame':	1428 obs. of  3 variables:
#   $ abundances:List of 1428
# ..$ : num [1:12, 1:178] 12 taxa, 178 size classes

Groups <- read_csv("../../NumericalModel/Multi-ZooSizeSpectrumModel/TestGroups.csv") # Load in functional group information
Groups <- pull(Groups,species)
# 
# Get Lat and Lon from Numerical Model. This will be the basis for the statistical model comparison
enviro_data <- read_csv("../../NumericalModel/Multi-ZooSizeSpectrumModel/enviro_5d_3.csv") # Load in environmental data
tt = which(enviro_data$sst < 1) # remove squares where sst < 1 (poor coverage)
enviro_data = enviro_data[-tt,]

df <- enviro_data %>%
  dplyr::select("x","y") %>%
  rename("Lon" = x, "Lat" = y)

# enviro_data <- enviro_data %>%
#   dplyr::select(-c( "x", "y"))

# Loop through all output and build a raster for all groups for all gridsquares
# Assign each spcies to a column, sum it, and store.
for (i in 3:9){
  out <- map(ssm$abundances, `[`, i,)
  assign(Groups[i],unlist(map(out,sum)))
  df <- cbind(df, sapply(Groups[i], function(x) eval(parse(text = x))))
}

# df <- tibble(OmniCopepods, CarnCopepods, Chaetognaths, Euphausiids, Larvaceans, Jellyfish, Salps)

num_brick <- rasterFromXYZ(df)  #Convert first two columns as lon-lat and third as value 
names(num_brick) <- Groups[3:9]
saveRDS(num_brick, "ModelOutput/GlobalLayers/Num_Brick.rds")

#######
#######
w=10^(seq(from=-10.7, to = 7, 0.1)) # Vector to convert abundance to weights

df <- enviro_data %>%
  dplyr::select("x","y") %>%
  rename("Lon" = x, "Lat" = y)

# Loop through all output and build a raster for all groups for all gridsquares
# Assign each spcies to a column, sum it, and store.
for (i in 3:9){
  out <- map(ssm$abundances, `[`, i,)
  assign(Groups[i],unlist(map(out, ~ sum(. * w))))
  df <- cbind(df, sapply(Groups[i], function(x) eval(parse(text = x))))
}


num_brick <- rasterFromXYZ(df)  #Convert first two columns as lon-lat and third as value 
names(num_brick) <- Groups[3:9]
saveRDS(num_brick, "ModelOutput/GlobalLayers/Num_Brick_Biomass.rds")


############# ############# ############# ############# ############# 
###### BUILD/LOAD IDENTICAL RASTER FROM THE STATISTICAL MODELS  #####
############# ############# ############# ############# ############# 
latlonCRS <- CRS("+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m")

for (i in 3:9){
  # Load global layer from statistical model and use calc() to apply functions over a raster object
  dat <- readRDS(paste0("ModelOutput/GlobalLayers/StatModel_Layer_",Groups[i],".rds"))
  mn <- calc(dat, fun = mean) # Don't remove NAs so that the polar cells are removed, rather than skewed towards summer
  names(mn) = Groups[i]
  if (i == 3) {stat_brick <- brick(resample(mn, num_brick, method='bilinear'))} else 
  {brick <- brick(resample(mn, num_brick, method='bilinear')) 
  stat_brick <- addLayer(stat_brick,brick)}
  rm(dat,mn,brick)
}
stat_brick <- brick(stat_brick)
saveRDS(stat_brick, "ModelOutput/GlobalLayers/Stat_Brick.rds")

# 
