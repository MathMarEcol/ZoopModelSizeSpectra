
library(tidyverse)
library(lubridate)
library(ggpubr)
library(ncdf4)
library(maptools)
library(raster)
library(rgdal)
data(wrld_simpl)

quick <- 1

EnviroDir <- "Data/EnvironmentalData"

########### Import CHLOROPHYLL AND SST ###########
SST_files <- list.files(path = EnviroDir, pattern = "SST_raster_", recursive = FALSE, full.names = TRUE)
Chl_files <- list.files(path = EnviroDir, pattern = "Chl_raster_", recursive = FALSE, full.names = TRUE)

# Ensure they are in order
SST_files <- sort(SST_files)
Chl_files <- sort(Chl_files)

SST_files <- SST_files[2:13]
Chl_files <- Chl_files[2:13]

lab <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
         "Sep", "Oct", "Nov", "Dec")

# Define themes for plots 
theme_opts <- list(theme(
  panel.grid.minor = element_blank(),
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

myplots <- list()
for(i in 1:12){
  
  rasr <- readRDS(SST_files[i])
  ras_df <- raster::as.data.frame(rasr,xy = T)
  
  col <- rev(rasterImage::colorPalette(n = 14, type = "jet.colors")) # using rasterImage pckge
  
  ## Do the plotting
  gg <- ggplot() +
    geom_tile(data = ras_df, aes(x = x, y = y, fill = SST)) +
    # geom_polygon(fill = "white") +
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group)) + 
    {if (quick==1) coord_quickmap() else coord_map(projection = "mollweide") } +
    scale_fill_gradientn(name =  expression(paste("SST (°C)")), colours = rev(col),
                         position = "bottom", na.value = "grey50",
                         limits = c(-2, 34)) +
    guides(fill = guide_colourbar(barwidth = 6, barheight = 0.5, 
                                  label.theme = element_text(size = 8),
                                  title.theme = element_text(size = 8))) +
    theme_opts +
    scale_alpha(range = c(-0, 0.5))
  
  myplots[[i]] <- gg
  rm(rasr, ras_df)
}

graphics.off()
ggarrange(plotlist=myplots, ncol = 3, nrow = 4, 
          labels=lab, hjust = -0.1,
          font.label = list(size = 10, color = "black", face = "bold", family = NULL))

ggsave("Data/EnvironmentalData/SST.png", device = "png", dpi=300, scale=1,
       width = 210, height = 297, units = "mm")


########## ########## ########## ##########
########## ########## ########## ##########
########## ########## ########## ##########

myplots <- list()

for(i in 1:12){
  
  # datmth <- datNA %>% 
  #   filter(Month==i)
  
  rasr <- readRDS(Chl_files[i])
  ras_df <- raster::as.data.frame(rasr,xy = T)
  
  ras_df <- ras_df %>% 
    mutate(Chl <- log10(Chl),
           Chl = replace(Chl,
                         Chl > 1,1))
  
  col <- rev(rasterImage::colorPalette(n = 14, type = "jet.colors")) # using rasterImage pckge
  
  ## Do the plotting
  gg <- ggplot() +
    geom_tile(data = ras_df, aes(x = x, y = y, fill = Chl)) +
    # geom_polygon(fill = "white") +
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group)) + 
    scale_fill_gradientn(name =  expression(paste("mg m3")), colours = rev(col),
                         position = "bottom", na.value = "grey50",
                         limits = c(-2, 1)) +
    guides(fill = guide_colourbar(barwidth = 6, barheight = 0.5, 
                                  label.theme = element_text(size = 8),
                                  title.theme = element_text(size = 8))) +
    {if (quick==1) coord_quickmap() else coord_map(projection = "mollweide") } +
    theme_opts +
    scale_alpha(range = c(-0, 0.5))
  # +
  #   geom_point(data = datmth, aes(x=Longitude, y = Latitude), color = "darkred", size = 1)
  # 
  myplots[[i]] <- gg
  rm(rasr, ras_df)
}


graphics.off()
ggarrange(plotlist=myplots, ncol = 3, nrow = 4, labels=lab, hjust = -0.1,
          font.label = list(size = 10, color = "black", face = "bold", family = NULL))

ggsave("Data/EnvironmentalData/Chl.png", device = "png", dpi=300, scale=1,
       width = 210, height = 297, units = "mm")


########## ########## ########## ##########
########## ########## ########## ##########
########## ########## ########## ##########
rasr <- readRDS("Data/EnvironmentalData/Bathy_raster_1Deg.rds")
ras_df <- raster::as.data.frame(rasr,xy = T)
# ras_df <- ras_df %>% 
#   mutate(Elevation.relative.to.sea.level = replace(Elevation.relative.to.sea.level,
#                                                    Elevation.relative.to.sea.level > 0,0),
#          Elevation.relative.to.sea.level = abs(Elevation.relative.to.sea.level))

col <- rev(rasterImage::colorPalette(n = 14, type = "jet.colors")) # using rasterImage pckge

## Do the plotting
gg <- ggplot() +
  geom_tile(data = ras_df, aes(x = x, y = y, fill = Bathy)) +
  # geom_polygon(fill = "white") +
  geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group)) + 
  scale_fill_gradientn(name =  expression(paste("m")), colours = rev(col),
                       position = "bottom", na.value = "grey50") +
  guides(fill = guide_colourbar(barwidth = 6, barheight = 0.5, 
                                label.theme = element_text(size = 8),
                                title.theme = element_text(size = 8))) +
  {if (quick==1) coord_quickmap() else coord_map(projection = "mollweide") } +
  theme_opts +
  scale_alpha(range = c(-0, 0.5))

ggsave("Data/EnvironmentalData/Bathy.png", device = "png", dpi=300, scale=1,
       height = 210, width = 297, units = "mm")
rm(rasr, ras_df)

##### Mean SST
rasr <- readRDS("Data/EnvironmentalData/SST_raster_00Mean_1Deg.rds")
ras_df <- raster::as.data.frame(rasr,xy = T)
col <- rev(rasterImage::colorPalette(n = 14, type = "jet.colors")) # using rasterImage pckge

## Do the plotting
gg <- ggplot() +
  geom_tile(data = ras_df, aes(x = x, y = y, fill = SST)) +
  # geom_polygon(fill = "white") +
  geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group)) + 
  {if (quick==1) coord_quickmap() else coord_map(projection = "mollweide") } +
  scale_fill_gradientn(name =  expression(paste("SST (°C)")), colours = rev(col),
                       position = "bottom", na.value = "grey50",
                       limits = c(-2, 34)) +
  guides(fill = guide_colourbar(barwidth = 6, barheight = 0.5, 
                                label.theme = element_text(size = 8),
                                title.theme = element_text(size = 8))) +
  theme_opts +
  scale_alpha(range = c(-0, 0.5))

ggsave("Data/EnvironmentalData/SST_Mn.png", device = "png", dpi=300, scale=1,
       height = 210, width = 297, units = "mm")
rm(rasr, ras_df)


##### Mean Chl
rasr <- readRDS("Data/EnvironmentalData/Chl_raster_00Mean_1Deg.rds")
ras_df <- raster::as.data.frame(rasr,xy = T)
ras_df <- ras_df %>% 
  mutate(Chl <- log10(Chl),
         Chl = replace(Chl,Chl > 1,1))

col <- rev(rasterImage::colorPalette(n = 14, type = "jet.colors")) # using rasterImage pckge

## Do the plotting
gg <- ggplot() +
  geom_tile(data = ras_df, aes(x = x, y = y, fill = Chl)) +
  # geom_polygon(fill = "white") +
  geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group)) + 
  {if (quick==1) coord_quickmap() else coord_map(projection = "mollweide") } +
  scale_fill_gradientn(name =  expression(paste("mg m3")), colours = rev(col),
                       position = "bottom", na.value = "grey50",
                       limits = c(-2, 1)) +
  guides(fill = guide_colourbar(barwidth = 6, barheight = 0.5, 
                                label.theme = element_text(size = 8),
                                title.theme = element_text(size = 8))) +
  theme_opts +
  scale_alpha(range = c(-0, 0.5))

ggsave("Data/EnvironmentalData/Chl_Mn.png", device = "png", dpi=300, scale=1,
       height = 210, width = 297, units = "mm")
rm(rasr, ras_df)