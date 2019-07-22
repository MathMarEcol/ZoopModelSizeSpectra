fPlotGlobalComparison <- function(layer){
  
  source("fBuildRasterBrick.R")
  source("fWorldBordersMoll.R")
  
  mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")
  CRS <- CRS("+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")
  
  # Define themes for plots 
  theme_opts <- list(theme(panel.grid.minor = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.background = element_blank(),
                           panel.border = element_blank(),
                           plot.background = element_rect(fill="white"),
                           plot.title = element_text(hjust = 0.5),
                           axis.line = element_blank(),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks = element_blank(),
                           axis.title.x = element_blank(),
                           axis.title.y = element_blank(),
                           legend.title = element_text(size=8),
                           legend.text = element_text(size=8),
                           legend.position= "bottom", 
                           # legend.justification = "centre",
                           ))
  
  # Load the World Borders in the Moll projection
  wb_df_moll <- fWorldBordersMoll()
  data_moll <- projectRaster(layer, crs = mollCRS, over = TRUE)
  data_moll_df <- as.data.frame(rasterToPoints(data_moll))
  
  ## Colormap options
  # col <- rainbow(14) # base R
  # col <- rev(rasterImage::colorPalette(n = 14, type = "bwr")) # using rasterImage pckge
  # col <- colorRampPalette(c("red", "yellow", "blue"))(14) # base R
  col <- rev(rasterImage::colorPalette(n = 20, type = "jet.colors")) # using rasterImage pckge
  
  # ## Do the plotting
  gg <- ggplot() +
    geom_tile(data = data_moll_df, aes(x = x, y = y, fill = layer)) +
    geom_polygon(fill = "white") +
    geom_polygon(data = wb_df_moll, aes(x = long, y = lat, group = group)) +
    scale_fill_gradientn(name =  expression(paste("ind. m"^-3)), colours = rev(col)) + #,
                         # limits = clim, breaks = cticks,
                         # labels = clabel, position = "bottom") +
    theme_opts +
    theme(plot.margin = unit(c(1,1,0,1), "cm")) +
    guides(fill = guide_colorbar(barheight = 0.3)) +
    scale_alpha(range = c(-0, 0.5)) +
    #coord_map(projection = "mollweide")
    coord_quickmap()
#barwidth
return(gg)  
}

