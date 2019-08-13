# Last Updated 29th May 2019

fWorldBorders <- function(){
  library(sf)
  shapef <- "~/Dropbox/EcosystemEfficiency/StatisticalModel/ZoopAbundanceModels/WorldBorders/TM_WORLD_BORDERS-0.3.shp"
  
  # Change projection of input shapefile to mollweide
  wb_sf <- read_sf(shapef)
  wb_sp <- as(wb_sf, "Spatial")
  CRS <- CRS("+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")
  # mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84") 
  # mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m")
  # 
  wb_moll <- spTransform(x = wb_sp, CRS)
  
  # Mollweide projection of World Borders
  wb_df <- fortify(wb_moll)  
  wb_df$group <- as.numeric(wb_df$group)
  wb_df$hole <- as.numeric(wb_df$hole)
  
  return(wb_df)
}