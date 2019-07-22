library(ggpubr)
library(sp)
library(tidyverse)
library(sf)

source("fWorldBordersMoll.R")

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

dat <- readRDS("LatestDatabaseOuput_Final_Enviro.rds")
taxa <- unique(dat$Group)

mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m")
latlonCRS <- CRS("+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m")

# Load the World Borders in the Moll projection
wb_df <- fWorldBordersMoll()

myplots <- list() # Initialise plot list

for (i in 1:length(taxa)) {
  sub <- dat %>% filter(Group==taxa[i]) %>%
    dplyr::select(c(Longitude, Latitude))
  
  sub_spdf <- SpatialPoints(coords = sub, proj4string = latlonCRS)
  sub_spdf = spTransform(sub_spdf,mollCRS)
  sub_df <- as.data.frame(sub_spdf)
  
  ##################################
  ## Plot data distribution
  fname <- paste0("Figures/DataDistribution/DataDistribution_",taxa[i],".png")
  gg <- 
    ggplot() + 
    geom_polygon(fill="white") + 
    geom_polygon(data = wb_df, aes(x = long, y = lat, group = group)) +
    geom_point(data = sub_df, aes(x = Longitude, y = Latitude), size=0.1, color="blue") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(taxa[i]) + 
    theme_opts +
    coord_map(projection = "mollweide")
  
    ggsave(gg,filename = fname, width = 297, height = 210, 
           units = "mm", dpi=300)
  
  myplots[[i]] <- gg
  
  rm(sub)
}


graphics.off()
ggarrange(plotlist=myplots, ncol = 3, nrow = 3)

ggsave("Figures/DataDistribution/DataDistribution_All.png", dpi=300,
       width = 210, height = 297, units = "mm")

