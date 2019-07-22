# Last Updated 29th May 2019

fWorldBordersMoll <- function(){
  library(sf)
  shapef <- "~/Dropbox/EcosystemEfficiency/StatisticalModel/ZoopAbundanceModels/WorldBorders/TM_WORLD_BORDERS-0.3.shp"
  
  # Change projection of input shapefile to mollweide
  wb_sf <- read_sf(shapef)
  wb_sp <- as(wb_sf, "Spatial")
  mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84") 
  mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m")
  
  wb_moll <- spTransform(x = wb_sp, mollCRS)
  
  # Mollweide projection of World Borders
  wb_df_moll <- fortify(wb_moll)  
  wb_df_moll$group <- as.numeric(wb_df_moll$group)
  wb_df_moll$hole <- as.numeric(wb_df_moll$hole)
  
  return(wb_df_moll)
}