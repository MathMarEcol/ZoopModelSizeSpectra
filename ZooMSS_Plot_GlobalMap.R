library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(raster)

source("fZooMSS_Xtras.R")

w <- 10^(seq(from = -10.7, to =  7, 0.1)) # Size bins of whole model

res <- read_rds("Output/res_20200429_Control.RDS")
groups <- read_csv("TestGroups.csv")
enviro <- read_rds("envirofull_20200317.RDS")

# Calculate biomass of each species for each cell
bioms <- fZooMSS_SpeciesBiomass(res,w)
bioms <- as.data.frame(t(as.data.frame(bioms)))
colnames(bioms) <- groups$Species # Add names
rownames(bioms) <- c() # Remove row names

bioms <- bind_cols(bioms, enviro) # Bind enviro data

# Which variable to plot
var = "Fish_Small"

# What to plot
df <- bioms %>%
  dplyr::select("Lon", "Lat", all_of(var))

# Now do the plotting
world <- ne_countries(scale = "medium", returnclass = "sf")
world_sf <- st_transform(world, crs = 54009) # Convert to Mollweide

df_raster <- rasterFromXYZ(df)  #Convert first two columns as lon-lat and third as value
df_poly <- rasterToPolygons(df_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
df_sf <- st_as_sf(df_poly, crs = 4326) # Convert to sf

st_crs(df_sf) <- 4326 # Set base projection
df_sf_moll <- st_transform(df_sf, crs = 54009) # Alter the CRS to mollweide

gg_world <- ggplot() +
  geom_sf(data = df_sf_moll, aes_string(fill = var), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey50") +
  scale_fill_gradient(low = "white",
                       high = "lightcoral",
                       space = "Lab",
                       na.value = "grey50",
                       aesthetics = "fill",
                       oob = scales::squish,
                       guide = guide_colourbar(title = paste0(var, " Biomass"),
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  theme(plot.background = element_rect(fill = NA),
        plot.margin = unit(c(0,0,0,0), "mm"),
        panel.grid.major = element_line(colour = "grey50", size = 0.2),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(20,"mm"),
        legend.key.height = unit(4,"mm"))


graphics.off()
x11(width = 8, height = 6)
gg_world
ggsave(paste0("Figures/Biomass_",var,".png"), dpi = 300)
