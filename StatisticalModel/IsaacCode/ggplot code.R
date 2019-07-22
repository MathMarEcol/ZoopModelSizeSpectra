

# data = list object. Each element of that list is a raster
# shapefile = world shapefile to plot (check the GCM Dropbox folder for the shapefile)
# outdir = directory to save the plot
# rcp = warming scenario
# yrs = time frame

# for palette i'm currently using different sources... you can try RColorBrewer package and define your style


rs.ratio2 <- function(data, shapefile, outdir, rcp, rrun = "r1i1p1", yrs) {
  
  library(raster)
  library(rasterImage)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Change projection of input shapefile to mollweide
    wb_sf <- read_sf(shapefile)
    wb_sp <- as(wb_sf, "Spatial")
    mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84") #"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    wb_moll <- spTransform(x = wb_sp, mollCRS)
  # Define themes for plots 
    theme_opts <- list(theme(panel.grid.minor = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.background = element_blank(),
                             plot.background = element_rect(fill="#e6e8ed"),
                             panel.border = element_blank(),
                             axis.line = element_blank(),
                             axis.text.x = element_blank(),
                             axis.text.y = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank(),
                             plot.title = element_text(size=22)))
  # Loop to get plots by depth
    for (i in 1:length(data)) {
      # Isolate SD layer and transform values into cube root 
      single <- data[[i]]$layer
    # Reclassify layer into categories for plotting
      rng_cv <- c(-70, -50, 1, 
                  -50, -10, 2, 
                  -10, -5, 3, 
                  -5, -3, 4, 
                  -3, -2, 5, 
                  -2, -1, 6,
                  -1, 0, 7,
                  0, 1, 8,
                  1, 2, 9,
                  2, 3, 10,
                  3, 5, 11,
                  5, 10, 12,
                  10, 50, 13,
                  50, 70, 14)
      mt.cv <- matrix(rng_cv, ncol = 3, byrow = TRUE)
      cv <- c(" < -50", "-50 - -10", "-10 - -5", "-5 - -3", "-3 - -2", "-2 - -1", "-1 - 0",
              "0 - 1", "1 - 2", "2 - 3", "3 - 5", "5 - 10", "10 - 50", " > 50")
      single2 <- reclassify(single, mt.cv, include.lowest = FALSE) # VoCC magnitude categories
      
    # General file's names
      name.fig <- paste(i, rcp, rrun, yrs, sep = "_")
    # Mollweide projection of World Borders
      wb_df_moll <- fortify(wb_moll)  
      wb_df_moll$group <- as.numeric(wb_df_moll$group)
      wb_df_moll$hole <- as.numeric(wb_df_moll$hole)
    # Plotting mollweide projection world border + VoCC SD
      name.mag <- paste(outdir, "VoCC.Ratio_Scenarios", name.fig, ".pdf", sep = "")
      mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")
      test_moll <- projectRaster(single2, crs = mollCRS, over = TRUE)
      testing <- as.data.frame(rasterToPoints(test_moll))
      col <- rainbow(14) # base R
      col2 <- rev(rasterImage::colorPalette(n = 14, type = "bwr")) # using rasterImage pckge
      col3 <- colorRampPalette(c("red", "yellow", "blue"))(14) # base R
      ggplot() + 
        geom_tile(data = testing, aes(x = testing$x, y = testing$y, fill = testing$layer)) +
        geom_polygon(fill="white") + 
        geom_polygon(data = wb_df_moll, aes(x = long, y = lat, group = group)) +
        scale_fill_gradientn(name="", 
                             colours = rev(col3), # col[min(cv_rcl[], na.rm = TRUE):max(cv_rcl[], na.rm = TRUE)]
                             limits = c(1, 14), breaks = seq(1,14,1), labels = cv, position = "bottom") +
        scale_alpha(range = c(-0, 0.5)) + 
        theme_opts + 
        theme(legend.position="bottom") + 
        guides(fill = guide_colorbar(barheight = 2, barwidth = 50)) +
        ggsave(filename = name.mag, width=20, height=12, dpi=300)
  } 
  
}