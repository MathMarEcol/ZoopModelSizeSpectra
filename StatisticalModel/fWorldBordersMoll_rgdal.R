
# Last Updated 20th June 2019

library(rgdal)
library(broom)

fWorldBordersMoll <- function(){
  shapef <- "WorldBorders/TM_WORLD_BORDERS-0.3.shp"
  
  # Change projection of input shapefile to mollweide
  wb_sf <- readOGR(shapef)
  
  # mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84") 
  mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m")
  wb_moll <- spTransform(x = wb_sf, mollCRS)
  
  wb_sp <- as(wb_moll, "Spatial")
  
  # Mollweide projection of World Borders
  wb_df_moll <- tidy(wb_moll)  
  wb_df_moll$group <- as.numeric(wb_df_moll$group)
  wb_df_moll$hole <- as.numeric(wb_df_moll$hole)
  
  return(wb_df_moll)
}