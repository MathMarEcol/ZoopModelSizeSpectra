library(ggpubr)
library(sp)
library(tidyverse)
library(sf)

source("WorldBordersMoll.R")

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



dat <- read_csv("Data/COPEPOD_Database/copepod__4000000-compilation.csv", 
                na = c("null","----","-----","-99.000","n/a","-"))


mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")
CRS <- CRS("+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")

# Load the World Borders in the Moll projection
# wb_df_moll <- WorldBordersMoll()

## This code is temp bnecause i can't plot as Mollwider. When I can delete, and go back to the function above.
shapef <- "WorldBorders/TM_WORLD_BORDERS-0.3.shp"
wb_sf <- read_sf(shapef)
wb_sp <- as(wb_sf, "Spatial")
wb <- spTransform(x = wb_sp, CRS)

# Mollweide projection of World Borders
wb_df <- fortify(wb)  
wb_df$group <- as.numeric(wb_df$group)
wb_df$hole <- as.numeric(wb_df$hole)

# wb_moll <- spTransform(wb_moll, CRS)

# sub <- dat %>% filter(Group==taxa[i] & TotAbundance>0) %>%
sub <- dat %>%
  dplyr::select(c(LATITUDE,LONGITDE))

sub_spdf <- SpatialPointsDataFrame(coords = sub, data = sub, proj4string = CRS)
sub_df <- as.data.frame(sub_spdf)

##################################
## Plot data distribution
# fname <- paste0("Figures/DataDistribution/DataDistribution_",taxa[i],"_NoZeros.png")
fname <- "Figures/COPEPOD_DataDistribution.png"

gg <- 
  ggplot() + 
  geom_polygon(fill="white") + 
  geom_polygon(data = wb_df, aes(x = long, y = lat, group = group)) +
  geom_point(data = sub_df, aes(x = LONGITDE, y = LATITUDE), size=0.1, color="blue") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_opts +
  ggsave(filename = fname, width=10, height=6, dpi=300)



dat <- dat %>% filter(LONGITDE > -180 & LONGITDE < -80 & LATITUDE < 20 & UNITS == "#/m3")

dat <- read_csv("Data/COPEPOD_Database/copepod__4000000-compilation.csv", 
                na = c("null","----","-----","-99.000","n/a","-"))

dat <- dat %>% filter(`DATASET-ID` == "EASTROPAC" & UNITS == "#/m3")






