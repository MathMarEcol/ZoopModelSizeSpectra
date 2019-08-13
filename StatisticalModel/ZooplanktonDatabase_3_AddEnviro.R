library(tidyverse)
library(ncdf4)
library(yaImpute)
library(lubridate)
library(raster)

source("func/fAddEnviro.R")

EnviroDir <- "Data/EnvironmentalData"
bathy_file <- paste0(EnviroDir,"/","GEBCO_2014_2D.nc")

########### IMPORT ZOOPLANKTON DATA ###########
dat <- readRDS("DatabaseOutput/LatestDatabaseOuput_Final.rds")

dat <- dat %>% mutate(
  Day = replace(Day,is.na(Day), 15)) # Make all missing days to the middle of the month
  
########### ADD ENVIRO DATA ###########
# Function which adds satellite SST, Chl and Bathy data to any data frame with "Latitude" and "Longitude" variables
dat <- fAddEnviro(dat, EnviroDir, bathy_file)


########### BATHY ###########
dat <- dat %>% mutate(
  Bathy = replace(Bathy, Bathy >= -20, -20), # 0s and NAs must be where it is shallow - recode them to -20 m
  Bathy = replace(Bathy, Bathy < -7000, -7000),
  Bathy = replace(Bathy, is.na(Bathy), -20),
  Bathy = abs(Bathy))  # Convert all bathymetry to positive

summary(dat$Bathy)

# col <- rev(rasterImage::colorPalette(n = 14, type = "jet.colors")) # using rasterImage pckge
# ggplot(data = dat, aes(x = Longitude, y = Latitude)) +
#   geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='lightgrey', colour = "black") +
#   geom_point(aes(color=Bathy)) + 
#   scale_fill_gradientn(colours = rev(col))


########### CHL ###########
summary(dat$Chl) # 89632 NAs

dat <- dat %>% filter(!is.na(Chl)) %>%  # Remove NAs
  mutate(Chl = replace(Chl, Chl > 10, 10)) # Set max Chl

summary(dat$Chl)

hist(dat$Chl)
hist(log10(dat$Chl))

## THIS STILL NEEDS TO BE FIXED
# # Some are obviously on land though, so lets find out what the problem is....
# rasChl_df <- rasChl_df %>% 
#   mutate(Chlorophyll.Concentration..OCI.Algorithm <- log10(Chlorophyll.Concentration..OCI.Algorithm),
#          Chlorophyll.Concentration..OCI.Algorithm = replace(Chlorophyll.Concentration..OCI.Algorithm,
#                                                             Chlorophyll.Concentration..OCI.Algorithm > 1,1))
# 
# col <- rev(rasterImage::colorPalette(n = 14, type = "jet.colors")) # using rasterImage pckge
# 
# ggplot() +
#   geom_tile(data = rasChl_df, aes(x = x, y = y, fill = Chlorophyll.Concentration..OCI.Algorithm)) +
#   geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group)) + 
#   coord_fixed(ratio = 1.) + 
#   scale_fill_gradientn(name =  expression(paste("mg m3")), colours = rev(col),
#                        position = "bottom", na.value = "grey50",
#                        limits = c(-2, 1)) +
#   guides(fill = guide_colourbar(barwidth = 6, barheight = 0.5, 
#                                 label.theme = element_text(size = 8),
#                                 title.theme = element_text(size = 8))) +
#   scale_alpha(range = c(-0, 0.5)) + 
#   geom_point(data = proje, aes(x = Longitude, y = Latitude), color="darkred", size = 1)


########### SST ###########
summary(dat$SST) # 20484 NAs

dat <- dat %>% filter(!is.na(SST)) # Remove NAs
summary(dat$SST)
hist(dat$SST)


# # OK the majority of missing data are in the polar reasons - makes sense
# ggplot(data = datNA, aes(x = Longitude, y = Latitude)) +
#   geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") +
#   geom_point(color="red")

# 
# jelly <- dat %>% filter(Group == "Jellyfish" & Type == "Net")
# hist(log10(jelly$TotAbundance))
# hist(jelly$Mesh[jelly$Mesh<501])
# 
# salps <- dat %>% filter(Group == "Salps" & Type == "Net")
# hist(log10(salps$TotAbundance))
# hist(salps$Mesh[salps$Mesh<501])

saveRDS(dat, file = "DatabaseOutput/LatestDatabaseOuput_Final_Enviro.rds")


