### Author: Ryan Heneghan
### Updated by Jason Everett
### Date: May 2019
### Last Updated: 4th July 2019
### This script creates figures comparing output from the size spectrum model, with
### output from the statistical models for different zooplankton groups.

library(raster)
library(tidyverse)
library(ggpubr)

source("func/fPlotGlobalComparison.R")
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
##### ##### ##### STATISTICAL MODEL ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### 

sdf <- as.data.frame(stat_brick)

sdf <- sdf %>% 
  map_df(~replace(., is.nan(.), NA)) %>% 
  drop_na() %>% 
  mutate( #Euphausiids = replace(Euphausiids, Euphausiids>1e1, 1e1),
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
  ungroup()

ggplot(dat = sdf4, aes(x = Chl, y = MnAbundance)) +
  geom_point() +
  facet_wrap(vars(Taxa), scales = "free")

sdf5 <- sdf4 %>% 
  group_by(Chl) %>% 
  summarise(TotAbundance = sum(MnAbundance))

sdf6 <- left_join(sdf4, sdf5, by = "Chl")

sdf6 <- sdf6 %>% 
  mutate(Proportion = MnAbundance/TotAbundance,
         Chl = as.numeric(as.character(Chl)))

Span <- 0.5
LoessOmni <- loess(Proportion  ~ Chl, data = sdf6[sdf6$Taxa == "OmniCopepods",], span = Span) 
LoessCarn <- loess(Proportion  ~ Chl, data = sdf6[sdf6$Taxa == "CarnCopepods",], span = Span) 
LoessEuph <- loess(Proportion  ~ Chl, data = sdf6[sdf6$Taxa == "Euphausiids",], span = Span) 
LoessLarv <- loess(Proportion  ~ Chl, data = sdf6[sdf6$Taxa == "Larvaceans",], span = Span) 
LoessSalp <- loess(Proportion  ~ Chl, data = sdf6[sdf6$Taxa == "Salps",], span = Span) 
LoessChaet <- loess(Proportion  ~ Chl, data = sdf6[sdf6$Taxa == "Chaetognaths",], span = Span) 
LoessJelly <- loess(Proportion  ~ Chl, data = sdf6[sdf6$Taxa == "Jellyfish",], span = Span) 
dat <- as.data.frame(cbind(unique(sdf6$Taxa), rbind(LoessCarn$fitted, LoessChaet$fitted, LoessEuph$fitted, 
                                                    LoessJelly$fitted, LoessLarv$fitted, LoessOmni$fitted, LoessSalp$fitted)))

dat2 <- gather(dat, key = "Chl", value = "Proportion", V2:V23) 
dat2$Chl <- sdf6$Chl
dat2$Proportion <- as.numeric(dat2$Proportion)

ggplot(dat2, aes(x = Chl, y = Proportion, fill = V1)) + 
  scale_x_continuous(limits = c(-1.5001, 0.5001), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.001,1.001), expand = c(0, 0)) +
  xlab("log10(Chlorophyll)") + labs(fill = "Taxa") +
  geom_area() + theme_bw(base_size = 16)
ggsave("Figures/StatModel_Abund_Prop.png", device = "png", width = 10, height = 6, dpi = 1200)

##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### #####  NUMERICAL MODEL  ##### ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### 

ndf <- as.data.frame(num_brick)

ndf <- ndf %>% 
  mutate(Chl = replace(Chl, Chl > 0.5, 0.5)) %>% 
         drop_na() # Setting Chl to max log10(0.5)

ndf2 <- gather(ndf, key = "Taxa", value = "Abundance", -c(SST, Chl))

ndf3 <- ndf2 %>% 
  filter(!is.na(Chl) & !is.na(SST)) %>% 
  mutate(Chl = round(Chl,1), # Round to 2 decimal places
         Chl = as.factor(Chl)) %>% 
  dplyr::group_by(Chl, Taxa)

ndf4 <- ndf3 %>% 
  summarise(MnAbundance = mean(Abundance)) %>% 
  ungroup()

ndf5 <- ndf4 %>% 
  group_by(Chl) %>% 
  summarise(TotAbundance = sum(MnAbundance))

ndf6 <- left_join(ndf4, ndf5, by = "Chl")

ndf6 <- ndf6 %>% 
  mutate(Proportion = MnAbundance/TotAbundance,
         Chl = as.numeric(as.character(Chl)))

ggplot(ndf6, aes(x=Chl, y=Proportion, fill=Taxa)) + 
  geom_area()

ggsave("Figures/NumModel_Abund_Prop.png", device = "png", dpi=400, scale=1, units = "mm")


##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### #####   NUMERICAL MODEL - BIOMASS   ##### ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### 

biomass_brick <- readRDS("ModelOutput/GlobalLayers/Num_Brick_Biomass.rds")

biomass_brick <- addLayer(biomass_brick, mnChl2, mnSST2)


bdf <- as.data.frame(biomass_brick)
bdf <- bdf %>% 
  mutate(Chl = replace(Chl, Chl > 0.5, 0.5)) # Setting Chl to max log10(0.5)

bdf2 <- gather(bdf, key = "Taxa", value = "Biomass", -c(SST, Chl))

bdf3 <- bdf2 %>% 
  filter(!is.na(Chl) & !is.na(SST)) %>% 
  mutate(Chl = round(Chl,1), # Round to 2 decimal places
         Chl = as.factor(Chl)) %>% 
  dplyr::group_by(Chl)

bdf4 <- bdf3 %>% 
  summarise(TotalBio = sum(Biomass,na.rm = T),
            Larvaceans = sum(Biomass[Taxa == "Larvaceans"],na.rm = T)/TotalBio,
            OmniCopepods = sum(Biomass[Taxa == "OmniCopepods"],na.rm = T)/TotalBio,
            CarnCopepods = sum(Biomass[Taxa == "CarnCopepods"],na.rm = T)/TotalBio,
            Euphausiids = sum(Biomass[Taxa == "Euphausiids"],na.rm = T)/TotalBio,
            Chaetognaths = sum(Biomass[Taxa == "Chaetognaths"],na.rm = T)/TotalBio,
            Salps = sum(Biomass[Taxa == "Salps"],na.rm = T)/TotalBio,
            Jellyfish = sum(Biomass[Taxa == "Jellyfish"],na.rm = T)/TotalBio) %>% 
  ungroup()

bdf5 <- gather(bdf4, key = "Taxa", value = "Proportion", -c(TotalBio, Chl)) %>% 
  mutate(Chl = as.numeric(as.character(Chl)))

ggplot(bdf5, aes(x=Chl, y=Proportion, fill=Taxa)) + 
  geom_area()

ggsave("Figures/NumModel_Biomass_Prop.png", device = "png", dpi=400, scale=1, units = "mm")


### STATISTICAL BIOMASS PROPORTIONS
MeanMass <- c(0.002506187, 0.00015813, 0.001581155, 0.997662705, 0.0629469, 1.990560912, 1.990560912)

# Geometric Mean
MeanMass <- c(5.01187E-05, 3.16228E-06, 1E-05, 0.011220185, 0.000398107, 0.014125375, 0.630957344)

sdf <- as.data.frame(stat_brick)

sdf$Larvaceans <- 10^sdf$Larvaceans * MeanMass[1]
sdf$OmniCopepods <- 10^sdf$OmniCopepods * MeanMass[2]
sdf$CarnCopepods <- 10^sdf$CarnCopepods * MeanMass[3]
sdf$Euphausiids <- 10^sdf$Euphausiids * MeanMass[4]
sdf$Chaetognaths <- 10^sdf$Chaetognaths * MeanMass[5]
sdf$Salps <- 10^sdf$Salps * MeanMass[6]
sdf$Jellyfish <- 10^sdf$Jellyfish * MeanMass[7]


sdf <- sdf %>% 
  map_df(~replace(., is.nan(.), NA)) %>% 
  drop_na() %>% 
  mutate( #Euphausiids = replace(Euphausiids, Euphausiids>1e1, 1e1),
    Chl = replace(Chl, Chl > 0.5, 0.5)) # Setting Chl to max log10(0.5)

sdf2 <- gather(sdf, key = "Taxa", value = "Biomass", -c(SST, Chl))


sdf3 <- sdf2 %>% 
  filter(!is.na(Chl) & !is.na(SST)) %>% 
  mutate(Chl = round(Chl,1), # Round to 2 decimal places
         Chl = as.factor(Chl)) %>% 
  dplyr::group_by(Chl, Taxa)

sdf4 <- sdf3 %>% 
  summarise(MnBiomass = mean(Biomass)) %>% 
  ungroup()

sdf5 <- sdf4 %>% 
  group_by(Chl) %>% 
  summarise(TotBiomass = sum(MnBiomass))

sdf6 <- left_join(sdf4, sdf5, by = "Chl")

sdf6 <- sdf6 %>% 
  mutate(Proportion = MnBiomass/TotBiomass,
         Chl = as.numeric(as.character(Chl)))

ggplot(sdf6, aes(x=Chl, y=Proportion, fill=Taxa)) + 
  geom_area()

ggsave("Figures/StatModel_Biomass_Prop.png", device = "png", dpi=400, scale=1, units = "mm")


###




# 
# ## Do Comparison Plot 1
# pl <- list()
# lab <- list()
# let <- c("A)", "B)", "C)", "D)", "E)","F)", 
#          "G)", "H)", "I)", "J)", "K)", "L)")
# cnt <- 1 # counter
# 
# for (i in 1:4){
#   ## STATISTICAL MODEL PLOTS
#   sb <- subset(stat_brick, paste0(Groups[i+2]), drop=TRUE)
#   names(sb) <- "layer" # Change so name is always the same
#   crs(sb) <- latlonCRS
#   pl[[cnt]] <- fPlotGlobalComparison(sb)
#   lab[[cnt]] <- paste0(let[cnt], " Stat: ",Groups[i+2])
#   cnt <- cnt + 1
#   
#   ## NUMERICAL MODEL PLOTS
#   nb <- subset(num_brick, paste0(Groups[i+2]), drop=TRUE)
#   saveRDS(nb, file = paste0("ModelOutput/GlobalLayers/NumModel_Layer_",Groups[i+2],".rds"))
#   
#   names(nb) <- "layer" # Change so name is always the same
#   crs(nb) <- latlonCRS
#   nb <- calc(nb, log10)  
#   pl[[cnt]] <- fPlotGlobalComparison(nb)
#   lab[[cnt]] <- paste0(let[cnt], " Num: ",Groups[i+2])
#   cnt <- cnt + 1
#   
#   ## CORRELATION MODEL PLOTS
#   sb_df <- as.data.frame(sb)
#   nb_df <- as.data.frame(nb)
#   dat <- cbind(nb_df,sb_df)
#   names(dat) <- c("Statistical", "Numerical")
#   
#   pl[[cnt]] <- ggplot(dat, aes(Statistical, Numerical)) + 
#     geom_point() + 
#     geom_smooth(method = "lm", se = TRUE) +
#     xlab("Statistical") + 
#     ylab("Theoretical") +
#     theme_bw() + 
#     theme(aspect.ratio=1,
#           axis.title = element_text(size=8),
#           axis.text = element_text(size=8))
#   lab[[cnt]] <- paste0(let[cnt], " Corr: ",Groups[i+2])
#   cnt <- cnt + 1
#   
#   rm(sb, nb, dat)
# }
# ggarrange(plotlist=pl, ncol = 3, nrow = 4, 
#           widths = c(1.8, 1.8, 1), labels=lab, hjust = -0.1,
#           font.label = list(size = 10, color = "black", face = "bold", family = NULL))
# 
# ggsave("Figures/Comparison1.png", device = "png", dpi=400, scale=1,
#        width = 210, height = 297, units = "mm")
# 
# 
# ## Do Plot 2
# pl <- list()
# lab <- list()
# let <- c("A)", "B)", "C)", "D)", "E)", 
#          "F)", "G)", "H)", "I)", "J)", "K)", "L)")
# cnt <- 1 # counter
# 
# for (i in 5:7){
#   ## STATISTICAL MODEL PLOTS
#   sb <- subset(stat_brick, paste0(Groups[i+2]), drop=TRUE)
#   names(sb) <- "layer" # Change so name is always the same
#   crs(sb) <- latlonCRS
#   pl[[cnt]] <- fPlotGlobalComparison(sb)
#   lab[[cnt]] <- paste0(let[cnt], " Stat: ",Groups[i+2])
#   cnt <- cnt + 1
#   
#   
#   ## NUMERICAL MODEL PLOTS
#   nb <- subset(num_brick, paste0(Groups[i+2]), drop=TRUE)
#   saveRDS(nb, file = paste0("ModelOutput/GlobalLayers/NumModel_Layer_",Groups[i+2],".rds"))
#   
#   names(nb) <- "layer" # Change so name is always the same
#   crs(nb) <- latlonCRS
#   nb <- calc(nb, log10)  
#   pl[[cnt]] <- fPlotGlobalComparison(nb)
#   lab[[cnt]] <- paste0(let[cnt], " Num: ",Groups[i+2])
#   cnt <- cnt + 1
#   
#   ## CORRELATION MODEL PLOTS
#   sb_df <- as.data.frame(sb)
#   nb_df <- as.data.frame(nb)
#   dat <- cbind(nb_df,sb_df)
#   names(dat) <- c("Statistical", "Numerical")
#   
#   pl[[cnt]] <- ggplot(dat, aes(Statistical, Numerical)) + 
#     geom_point() + 
#     geom_smooth(method = "lm", se = TRUE) +
#     xlab("Statistical") + 
#     ylab("Theoretical") +
#     theme_bw() + 
#     theme(aspect.ratio=1,
#           axis.title = element_text(size=8),
#           axis.text = element_text(size=8))
#   lab[[cnt]] <- paste0(let[cnt], " Corr: ",Groups[i+2])
#   cnt <- cnt + 1
#   
#   rm(sb, nb, dat)
# }
# 
# ggarrange(plotlist=pl, ncol = 3, nrow = 4, 
#           widths = c(1.8, 1.8, 1), labels=lab, hjust = -0.1,
#           font.label = list(size = 10, color = "black", face = "bold", family = NULL))
# 
# ggsave("Figures/Comparison2.png", device = "png", dpi=400, scale=1,
#        width = 210, height = 297, units = "mm")
