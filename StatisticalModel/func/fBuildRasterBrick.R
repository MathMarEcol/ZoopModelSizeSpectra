fBuildRasterBrick <- function(NthMonth, maxm, res){
  
  library(lubridate)
  
  NthDOY <- lubridate::yday(parse_date_time(paste0("15 ", gsub("\\d", "", NthMonth), " 2010") , "dmy")) # Convert month to DOY assuming year 2010
  
  if (res == "one"){
    SST <- readRDS(paste0("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZoopAbundanceModels/Data/EnvironmentalData/SST_raster_",NthMonth,"_1Deg.rds"))
    }
  if(res =="onetenth"){
    SST <- readRDS(paste0("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZoopAbundanceModels/Data/EnvironmentalData/SST_raster_",NthMonth,"_oneTenthDeg.rds"))
  }
  if (sum(names(maxm)=="SST") > 0) { # Replace large values we didn't model
    SST <- reclassify(SST, c(maxm$SST,Inf,maxm$SST))
  }
  
  if (res == "one"){
    Chl <- readRDS(paste0("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZoopAbundanceModels/Data/EnvironmentalData/Chl_raster_",NthMonth,"_1Deg.rds"))
  }
  if(res =="onetenth"){
  Chl <- readRDS(paste0("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZoopAbundanceModels/Data/EnvironmentalData/Chl_raster_",NthMonth,"_oneTenthDeg.rds"))
  }
  if (sum(names(maxm)=="Chl") > 0) { # Replace large values we didn't model
    Chl <- reclassify(Chl, c(maxm$Chl,Inf,maxm$Chl))
  }
  
  
  if (res == "one"){
    Bathy <- readRDS("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZoopAbundanceModels/Data/EnvironmentalData/Bathy_raster_1Deg.rds")
  }
  if(res =="onetenth"){
    Bathy <- readRDS("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZoopAbundanceModels/Data/EnvironmentalData/Bathy_raster_oneTenthDeg.rds")
  }
if (sum(names(maxm)=="Bathy") > 0) { # Replace large values we didn't model
    Bathy <- reclassify(Bathy, c(maxm$Bathy,Inf,maxm$Bathy))
    Bathy <- reclassify(Bathy, c(-Inf,20,20))
  }
  
  # Create DOY2 Lookup table for Nth v Sth Hemi
  DOY_lookUp <- matrix(data = 1:366, nrow = 366, ncol = 2) 
  DOY_lookUp[,2] = c(184:366, 1:183)
  
  # Create a DOY2 raster layer
  DOY2 <- Bathy # Bathy already set up as raster
  names(DOY2) <- "DOY2"
  
  DOY2[1:90,] <- NthDOY
  DOY2[91:180,] <- DOY_lookUp[NthDOY,2]
  
  HarmDOY <- (DOY2/365)*2*pi
  
  # Prediction model
  # Make a raster brick with the 3-remotely-sensed datasets in there
  data_brick <- brick(Chl, SST, Bathy, HarmDOY)
  names(data_brick) <- c("Chl", "SST", "Bathy", "HarmDOY") # Rename raster brick layers so same as in lm
  
  return(data_brick)
  
}