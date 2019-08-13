fPlotGlobalGAMM <- function(NthMonth,taxa,df,excl,tframe, maxm){
  
  quick <- 1
  res <- "one" # Run at One Degree
  
  source("fBuildRasterBrick.R")
  source("fWorldBordersMoll.R")
  
  mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")
  CRS <- CRS("+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")
  
  myplots <- list() # Initialise plot list
  
  if (tframe=="month"){
    if (taxa == "OmniCopepods"){ticks <- c(100, 300, 1000, 3000, 10000)}
    if (taxa == "CarnCopepods"){ticks <- c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)}
    if (taxa == "Jellyfish"){ticks <- c(0.00003, 0.0001, 0.0003, 0.001)}
    if (taxa == "Chaetognaths"){ticks <- c(0.003, 0.01, 0.03, 0.1, 0.3, 1)}
    if (taxa == "Salps"){ticks <- c(0.001, 0.003, 0.01, 0.03)}
    if (taxa == "Larvaceans"){ticks <- c(0.1, 0.3, 1, 3, 10)}
    if (taxa == "Euphausiids"){ticks <- c(0.003, 0.01, 0.03, 0.1)}
  }
  
  if (tframe=="annual"){
    if (taxa == "OmniCopepods"){ticks <- c(5, 10, 40)}
    if (taxa == "CarnCopepods"){ticks <- c(0.1, 0.3, 1, 3, 10, 30)}
    if (taxa == "Euphausiids"){ticks <- c(0.1, 0.3, 1)}
    if (taxa == "Chaetognaths"){ticks <- c(0.7, 1.5, 3, 10, 20)}
    if (taxa == "Jellyfish"){ticks <- c(0.1, 0.2, 0.4, 1, 3, 10)}
    if (taxa == "Salps"){ticks <- c(0.03, 0.1, 0.3, 0.5, 1)}
    if (taxa == "Larvaceans"){ticks <- c(0.3, 1, 3, 10, 40)}
  }
  
  cticks <- log10(ticks)
  clabel <- as.character(ticks)
  clim <- c(cticks[1], last(cticks))
  
  # Define themes for plots 
  theme_opts <- list(theme(panel.grid.minor = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.background = element_blank(),
                           plot.background = element_rect(fill="white"),
                           plot.title = element_text(hjust = 0.5, size = rel(0.5)),
                           panel.border = element_blank(),
                           axis.line = element_blank(),
                           axis.text.x = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks = element_blank(),
                           axis.title.x = element_blank(),
                           axis.title.y = element_blank(),
                           legend.position = "bottom",
                           legend.text = element_text(size = rel(0.5))))
  
  for (i in 1:length(NthMonth)) {
    fname <- paste0("Figures/GlobalMaps/GlobalModel_Mth",sprintf('%02g', i),"_",taxa,".png")
    
    # The data brick of satellite data, and DOY
    data_brick <- fBuildRasterBrick(NthMonth[i], maxm, res)
    
    if (taxa == "OmniCopepods"){
      saveRDS(data_brick,paste0("ModelOutput/GlobalLayers/Enviro_Data_Brick_",NthMonth[i],".rds"))
      }
    
    # mdl <- readRDS(paste0("ModelOutput/gamm_",taxa,".rds"))
    mdl <- readRDS(paste0("ModelOutput/lmer_",taxa,"_log.rds"))
    
    # layer <- predict(data_brick, mdl, 
    #                  const = df, exclude = excl)

    layer <- predict(data_brick, mdl, 
                     const = df, re.form=NA)
    
        
    names(layer) <- paste0(gsub("\\d", "", NthMonth[i]))
    
    if (i == 1){
      stat_brick <- brick(layer)} else{
        stat_brick <- addLayer(stat_brick,layer)
      }
  }
  saveRDS(stat_brick, file = paste0("ModelOutput/GlobalLayers/StatModel_Layer_",taxa,".rds"))
  
  # Load the World Borders in the Moll projection
  wb_df_moll <- fWorldBordersMoll()
  
  if (tframe == "annual") {
    stat_brick  <- calc(stat_brick, fun = mean)
    n <- 1
  }
  if (tframe=="month"){
    n <- 12
  }
  
  # Convert data to Moll
  data_moll <- projectRaster(stat_brick, crs = mollCRS, over = TRUE)
  
  for (i in 1:n){
    # Convert data to dataframe
    data_moll_df <- as.data.frame(rasterToPoints(data_moll))
    
    ## Colormap options
    # col <- rainbow(14) # base R
    # col <- rev(rasterImage::colorPalette(n = 14, type = "bwr")) # using rasterImage pckge
    # col <- colorRampPalette(c("red", "yellow", "blue"))(14) # base R
    col <- rev(rasterImage::colorPalette(n = 14, type = "jet.colors")) # using rasterImage pckge
    
    # ## Do the plotting
    gg <- ggplot() +
      geom_tile(data = data_moll_df, aes_string(x = "x", y = "y", fill = names(data_moll_df[i+2]))) +
      geom_polygon(fill = "white") +
      geom_polygon(data = wb_df_moll, aes(x = long, y = lat, group = group)) +
      scale_fill_gradientn(name =  expression(paste("ind. m"^-3)), colours = rev(col),
                           limits = clim, breaks = cticks,
                           labels = clabel, position = "bottom",
                           na.value = "grey50") +
      {if (quick==1) coord_quickmap() else coord_map(projection = "mollweide") } +
      theme_opts +
      scale_alpha(range = c(-0, 0.5))
    
    if (tframe == "month"){
      gg <- gg +  ggtitle(paste0(gsub("\\d", "", NthMonth[i])," - ",taxa)) + 
        guides(fill = guide_colorbar(barheight = 0.5, barwidth = 10))
    }
    
    if (tframe == "annual"){
      gg <- gg +  ggtitle(taxa) + 
        guides(fill = guide_colorbar(barheight = 2, barwidth = 20)) + 
        theme(legend.text = element_text(size = rel(2)),
              legend.title = element_text(size = rel(2)),
              plot.title = element_text(size = rel(2)))
    }
    
    myplots[[i]] <- gg
  }
  
  
  if (tframe=="month"){
    graphics.off()
    x11(width = 15, height = 15)
    ggarrange(plotlist=myplots, ncol = 3, nrow = 4)
    ggsave(paste0("Figures/GlobalMaps/GlobalModel_All_",taxa,".png"), dpi=300, scale=1,
           width = 210, height = 297, units = "mm")
  } 
  
  if (tframe=="annual"){
    graphics.off()
    x11(width = 15, height = 15)
    ggsave(paste0("Figures/GlobalMaps/GlobalModel_Annual_",taxa,".png"), dpi=300, scale=1,
           width = 210, height = 297, units = "mm")
  } 
  
  return(myplots)
  
}
