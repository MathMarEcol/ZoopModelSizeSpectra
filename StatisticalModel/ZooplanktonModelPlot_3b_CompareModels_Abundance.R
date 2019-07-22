### Author: Ryan Heneghan
### Updated by Jason Everett
### Date: May 2019
### Last Updated: 4th July 2019
### This script creates figures comparing output from the size spectrum model, with
### output from the statistical models for different zooplankton groups.

library(raster)
library(tidyverse)
library(ggpubr)

source("fPlotGlobalComparison.R")
latlonCRS <- CRS("+proj=longlat +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m")

Groups <- read_csv("../../NumericalModel/Multi-ZooSizeSpectrumModel/TestGroups.csv") # Load in functional group information
Groups <- pull(Groups,species)

num_brick <- readRDS("ModelOutput/GlobalLayers/Num_Brick.rds")
stat_brick <- readRDS("ModelOutput/GlobalLayers/Stat_Brick.rds")

enviro_files <- sort(list.files(path = "ModelOutput/GlobalLayers", 
                                pattern = "Enviro_Data_Brick", 
                                recursive = FALSE, full.names = TRUE))
enviro_brick <- readRDS(paste0(enviro_files[1]))

Chl <- enviro_brick$Chl
SST <- enviro_brick$SST

for (i in 2:12){
  enviro_brick <- readRDS(paste0(enviro_files[i]))
  
  Chl <- addLayer(Chl,enviro_brick$Chl)
  SST <- addLayer(SST,enviro_brick$SST)
}

logChl = calc(Chl, log10)
mnChl = calc(logChl, mean)
names(mnChl) <- "Chl"

mnSST = calc(SST, mean)
names(mnSST) <- "SST"

mnChl2 <- resample(mnChl, num_brick, method='bilinear')
mnSST2 <- resample(mnSST, num_brick, method='bilinear')

num_brick <- addLayer(num_brick, mnChl2, mnSST2)
stat_brick <- addLayer(stat_brick, mnChl2, mnSST2)



##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### #####  NUMERICAL MODEL - ADUNDANCE  ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### 

ndf <- as.data.frame(num_brick)

ndf <- ndf %>% 
  map_df(~replace(., is.nan(.), NA)) %>% 
  drop_na() %>% 
  mutate(#Euphausiids = replace(Euphausiids, Euphausiids>1e1, 1e1),
    Chl = replace(Chl, Chl > 0.5, 0.5)) # Setting Chl to max log10(0.5)

ndf2 <- gather(ndf, key = "Taxa", value = "Abundance", -c(SST, Chl))

ndf3 <- ndf2 %>% 
  filter(!is.na(Chl) & !is.na(SST)) %>% 
  mutate(Chl = round(Chl,1), # Round to 2 decimal places
         Chl = as.factor(Chl)) %>% 
  dplyr::group_by(Chl, Taxa)

ndf4 <- ndf3 %>% 
  summarise(MnAbundance = mean(Abundance)) %>% 
  ungroup() %>% 
  mutate(Chl = as.numeric(as.character(Chl)))

ggplot(ndf4, aes(x=Chl, y=MnAbundance, fill=Taxa)) + 
  geom_area() +
  facet_wrap(vars(Taxa), ncol = 2, scales = "free")

ggsave("Figures/NumModel_Abund_RawAbund.png", 
       device = "png", dpi=400, scale=1, units = "mm")



##### ##### ##### ##### ##### ##### ##### ##### #####
##### ##### STATISTICAL MODEL - ADUNDANCE ##### #####
##### ##### ##### ##### ##### ##### ##### ##### #####

sdf <- as.data.frame(stat_brick)

sdf <- sdf %>%
  map_df(~replace(., is.nan(.), NA)) %>%
  drop_na() %>%
  mutate(#Euphausiids = replace(Euphausiids, Euphausiids>1e1, 1e1),
    Chl = replace(Chl, Chl > 0.5, 0.5)) # Setting Chl to max log10(0.5)

sdf2 <- gather(sdf, key = "Taxa", value = "Abundance", -c(SST, Chl))
sdf2 <- sdf2 %>%
  mutate(Abundance = 10^Abundance)

sdf3 <- sdf2 %>%
  filter(!is.na(Chl) & !is.na(SST)) %>%
  mutate(Chl = round(Chl,1), # Round to 2 decimal places
         Chl = as.factor(Chl)) %>%
  dplyr::group_by(Chl, Taxa)

sdf4 <- sdf3 %>%
  summarise(MnAbundance = mean(Abundance)) %>%
  ungroup() %>%
  mutate(Chl = as.numeric(as.character(Chl)))

ggplot(sdf4, aes(x=Chl, y=MnAbundance, fill=Taxa)) +
  geom_area() +
  facet_wrap(vars(Taxa), ncol = 2, scales = "free")

ggsave("Figures/StatModel_Abund_RawAbund.png", 
       device = "png", dpi=400, scale=1, units = "mm")




###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######

## Do Comparison Plot 1
pl <- list()
lab <- list()
let <- c("A)", "B)", "C)", "D)", "E)","F)",
         "G)", "H)", "I)", "J)", "K)", "L)")
cnt <- 1 # counter

for (i in 1:4){
  ## STATISTICAL MODEL PLOTS
  sb <- subset(stat_brick, paste0(Groups[i+2]), drop=TRUE)
  names(sb) <- "layer" # Change so name is always the same
  crs(sb) <- latlonCRS
  pl[[cnt]] <- fPlotGlobalComparison(sb)
  lab[[cnt]] <- paste0(let[cnt], " Stat: ",Groups[i+2])
  cnt <- cnt + 1

  ## NUMERICAL MODEL PLOTS
  nb <- subset(num_brick, paste0(Groups[i+2]), drop=TRUE)
  saveRDS(nb, file = paste0("ModelOutput/GlobalLayers/NumModel_Layer_",Groups[i+2],".rds"))

  names(nb) <- "layer" # Change so name is always the same
  crs(nb) <- latlonCRS
  nb <- calc(nb, log10)
  pl[[cnt]] <- fPlotGlobalComparison(nb)
  lab[[cnt]] <- paste0(let[cnt], " Num: ",Groups[i+2])
  cnt <- cnt + 1

  ## CORRELATION MODEL PLOTS
  sb_df <- as.data.frame(sb)
  nb_df <- as.data.frame(nb)
  dat <- cbind(nb_df,sb_df)
  names(dat) <- c("Statistical", "Numerical")

  pl[[cnt]] <- ggplot(dat, aes(Statistical, Numerical)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    xlab("Statistical") +
    ylab("Theoretical") +
    theme_bw() +
    theme(aspect.ratio=1,
          axis.title = element_text(size=8),
          axis.text = element_text(size=8))
  lab[[cnt]] <- paste0(let[cnt], " Corr: ",Groups[i+2])
  cnt <- cnt + 1

  rm(sb, nb, dat)
}
ggarrange(plotlist=pl, ncol = 3, nrow = 4,
          widths = c(1.8, 1.8, 1), labels=lab, hjust = -0.1,
          font.label = list(size = 10, color = "black", face = "bold", family = NULL))

ggsave("Figures/Comparison1.png", device = "png", dpi=400, scale=1,
       width = 210, height = 297, units = "mm")



###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### 

## Do Plot 2
pl <- list()
lab <- list()
let <- c("A)", "B)", "C)", "D)", "E)", 
         "F)", "G)", "H)", "I)", "J)", "K)", "L)")
cnt <- 1 # counter

for (i in 5:7){
  ## STATISTICAL MODEL PLOTS
  sb <- subset(stat_brick, paste0(Groups[i+2]), drop=TRUE)
  names(sb) <- "layer" # Change so name is always the same
  crs(sb) <- latlonCRS
  pl[[cnt]] <- fPlotGlobalComparison(sb)
  lab[[cnt]] <- paste0(let[cnt], " Stat: ",Groups[i+2])
  cnt <- cnt + 1
  
  
  ## NUMERICAL MODEL PLOTS
  nb <- subset(num_brick, paste0(Groups[i+2]), drop=TRUE)
  saveRDS(nb, file = paste0("ModelOutput/GlobalLayers/NumModel_Layer_",Groups[i+2],".rds"))
  
  names(nb) <- "layer" # Change so name is always the same
  crs(nb) <- latlonCRS
  nb <- calc(nb, log10)  
  pl[[cnt]] <- fPlotGlobalComparison(nb)
  lab[[cnt]] <- paste0(let[cnt], " Num: ",Groups[i+2])
  cnt <- cnt + 1
  
  ## CORRELATION MODEL PLOTS
  sb_df <- as.data.frame(sb)
  nb_df <- as.data.frame(nb)
  dat <- cbind(nb_df,sb_df)
  names(dat) <- c("Statistical", "Numerical")
  
  pl[[cnt]] <- ggplot(dat, aes(Statistical, Numerical)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = TRUE) +
    xlab("Statistical") + 
    ylab("Theoretical") +
    theme_bw() + 
    theme(aspect.ratio=1,
          axis.title = element_text(size=8),
          axis.text = element_text(size=8))
  lab[[cnt]] <- paste0(let[cnt], " Corr: ",Groups[i+2])
  cnt <- cnt + 1
  
  rm(sb, nb, dat)
}

ggarrange(plotlist=pl, ncol = 3, nrow = 4, 
          widths = c(1.8, 1.8, 1), labels=lab, hjust = -0.1,
          font.label = list(size = 10, color = "black", face = "bold", family = NULL))

ggsave("Figures/Comparison2.png", device = "png", dpi=400, scale=1,
       width = 210, height = 297, units = "mm")
