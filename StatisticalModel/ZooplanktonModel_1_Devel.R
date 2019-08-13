# ZooplanktonModelDevel.R
# Ant and Jase
# Created 11 April 2019
# Last Edited: 27 June 2019

## Memory Errors on Mac
# If you run into Memory errors do the following
# from Stack (https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
# ```cd ~
#    touch .Renviron
#    open .Renviron
#    Paste the following as the first line: R_MAX_VSIZE=100Gb 
#    Save
#    Exit ```

## NOTES
# 1. https://stats.stackexchange.com/questions/197952/two-methods-of-adding-random-effects-to-a-gam-give-very-different-results-why-i
# 2. Interpreting Q-Q plots. https://stats.stackexchange.com/questions/101274/how-to-interpret-a-qq-plot
# Our plots of heavy tailed on left and right (too many low values and too many high values). 
# Best are copepods because few zeros...

############################  Preliminaries ############################  
library(tidyverse)
library(visreg)
library(splines)
library(lme4)
library(MuMIn) # r2 for mixed models
library(lmerTest) # p-levels for mixed models # install.packages("lmerTest")

source("func/fHarmonic.R") # Harmonic function
source("func/fAddWeighting.R") # Because CPR is 90% of the data, use weights for individual Groups to give CPR 50% weight and Nets 50% weight
source("func/fPlotLMERs_Abundance.R") # Plot the Linear Mixed Effects Models
source("func/fPlotAbundanceLM.R") # Plot the Linear Mixed Effects Models
source("func/fPlotAbundanceLM_GearMesh.R") # Plot the Linear Mixed Effects Models
source("func/fSummariseLMERs.R") # Summaries calculating r2 and diagnostic plots

dat <- readRDS("DatabaseOutput/LatestDatabaseOuput_Final_Enviro.rds") # n = 1,009,208

dat <- dat %>% 
  mutate(
    HarmHour = (TimeLocal/24)*2*pi, # Convert to radians
    HarmDOY = (DOY2/365)*2*pi, # Convert to radians
    Latitude2 = abs(Latitude),
    TotAbundance = replace(TotAbundance, TotAbundance > 10000, 10000),
    Bathy = replace(Bathy, Bathy > 6000, 6000),
    Mid_Z = replace(Mid_Z, Mid_Z > 1000, 1000),
    Mesh = replace(Mesh, Mesh > 500, 500), 
    SST = replace(SST, SST > 30, 30),
    Chl = replace(Chl, Chl > 10, 10) # This has already been set to 10, but doing it again here for completeness and traceability
  ) %>% 
  mutate(
    Tow = case_when(Project == "SO-CPR" ~ "S1",
                    Project == "IMOS-CPR" ~ "S2",
                    Project == "SAHFOS" ~ "S3",
                    Tow == "H" ~ "H",
                    Tow == "V" ~ "V",
                    Tow == "O" ~ "O")) %>%
  dplyr::select(-Gear_Mesh) %>% unite(Gear_Mesh, Gear, Mesh, Tow, remove = FALSE) %>% # Remove old Gear_Mesh and calculate new Gear_Mesh
  mutate(Gear_Mesh = replace(Gear_Mesh, Gear_Mesh == "330_202.000007033348_O", "330_202_O")) %>% # Dodgy value...
  droplevels()

############################  GAMs ############################  
# RANDOM EFFECTS - Use best model here and add Longhurst and Gear. Maybe ShipCruise
# Gear as a random effect would be better as we don't need to choose a level in the maps
# HOW does the best model here (G1), work with the other functional groups

#### Define min_val for GAM for each Project within each Group 
Min <- dat %>% 
  filter(TotAbundance > 0) %>% 
  group_by(Group, Project) %>% 
  # summarise(Min = min(TotAbundance)) 
  summarise(Min = min(TotAbundance) / 2) # Use 1/2 lowest value

# Which samples that have zeros
Zeros <- dat %>%
  filter(TotAbundance == 0) %>% 
  group_by(Group, Project) %>% 
  summarise(N = n()) # Num of 0s

# Join all of the samples with those that have zeros
# Where there are NAs (i.e. no 0s) replace with 0 and then set Min to 0
# i.e. Those with no 0s, set log(Y+0)
All <- left_join(Min, Zeros, by = c("Group", "Project")) %>% 
  replace_na(list(N = 0)) %>% 
  mutate(Min = replace(Min, N == 0, 0))
dat <- left_join(dat, All, by = c("Group", "Project"))
rm(All, Min, Zeros)

# Add weights to give CPR and Nets equal weight in the analysis
dat_Omni <- fAddWeighting(dat, "OmniCopepods") %>% droplevels()
dat_Carn <- fAddWeighting(dat, "CarnCopepods")
dat_Chaet <- fAddWeighting(dat, "Chaetognaths")
dat_Euph <- fAddWeighting(dat, "Euphausiids")
dat_Jelly <- fAddWeighting(dat, "Jellyfish")
dat_Larv <- fAddWeighting(dat, "Larvaceans") %>% droplevels()
dat_Salp <- fAddWeighting(dat, "Salps")

########## Omnivorous Copepods ##########
m1_omni <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                  (1|Gear_Mesh),
                data = dat_Omni, weights = WtVec)
summary(m1_omni)
anova(m1_omni)
plot(m1_omni)
fSummariseLMERs(m1_omni, "Omnivores") # r2 = 77.7% (fixed = 3.9%)
fPlotLMERs_Abundance(m1_omni, "OmnivorousCopepods", "log")
saveRDS(m1_omni, file = "ModelOutput/lmer_OmniCopepods_log.rds")
rm(m1_omni)
graphics.off()

########## CarnCopepods ##########
m1_carn <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                  (1|Gear_Mesh),
                data = dat_Carn, weights = WtVec)
summary(m1_carn)
anova(m1_carn) # All significant
plot(m1_carn)
fSummariseLMERs(m1_carn, "Carnivores") # r2 = 75.3% (fixed = 7.6%)
fPlotLMERs_Abundance(m1_carn, "CarnivorousCopepods", "log")
saveRDS(m1_carn, file = "ModelOutput/lmer_CarnCopepods_log.rds")
rm(m1_carn)
graphics.off()

########## Chaetognaths ##########
m1_chaet <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Gear_Mesh),
           data = dat_Chaet, weights = WtVec)
summary(m1_chaet)
anova(m1_chaet) # All significant
plot(m1_chaet)
fSummariseLMERs(m1_chaet, "Chaetognaths") #r2 = 68.5% (fixed = 5.5%)
fPlotLMERs_Abundance(m1_chaet, "Chaetognaths", "log") 
saveRDS(m1_chaet, file = "ModelOutput/lmer_Chaetognaths_log.rds")
rm(m1_chaet)
graphics.off()

########## Larvaceans ##########
m1_larv <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Gear_Mesh),
           data = dat_Larv, weights = WtVec)
summary(m1_larv)
anova(m1_larv) # n.s. exp(-Mid_Z/1000):fHarmonic(HarmHour, k = 1)
fSummariseLMERs(m1_larv, "Larvaceans") # r2 = 65.9% (fixed = 3.9%)
fPlotLMERs_Abundance(m1_larv, "Larvaceans", "log")
saveRDS(m1_larv, file = "ModelOutput/lmer_Larvaceans_log.rds")
rm(m1_larv)
graphics.off()

########## Euphausiids ##########
# S1 = SO-CPR, S2=IMOS-CPR, S3 = SAHFOS
dat_Euph <-  dat_Euph %>% filter(Tow != "S3" & Tow != "S2" & Tow != "S1")
# dat_Euph <- fAddWeighting(dat_Euph, "Euphausiids")

# m1_euph <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
#              ns(Bathy, 3) + log10(Chl) + 
#              exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
#              (1|Gear_Mesh),
#            data = dat_Euph, weights = WtVec)
m1_euph <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                  (1|Gear_Mesh),
                data = dat_Euph)
summary(m1_euph)
anova(m1_euph) # All sig.
plot(m1_euph)
fSummariseLMERs(m1_euph, "Euphausiids") # r2 = 60.5% (fixed = 6.5%)
fPlotLMERs_Abundance(m1_euph, "Euphausiids", "log")
saveRDS(m1_euph, file = "ModelOutput/lmer_Euphausiids_log.rds")
rm(m1_euph)
graphics.off()

########## Salps ########## 
dat_Salp <-  dat_Salp %>% filter(Tow != "S3" & Tow != "S2" & Tow != "S1")
# m1_salp <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
#              ns(Bathy, 3) + log10(Chl) + 
#              exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
#              (1|Gear_Mesh),
#            data = dat_Salp, weights = WtVec)
m1_salp <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                  (1|Gear_Mesh),
                data = dat_Salp)
summary(m1_salp)
anova(m1_salp) # Mesh now p =0.09, but still +ve OLD: Mesh n.s. (p = 0.06), but Mesh positive
fSummariseLMERs(m1_salp, "Salps") # r2 = 54.2% (fixed = 6.8%)
fPlotLMERs_Abundance(m1_salp, "Salps", "log")
saveRDS(m1_salp, file = "ModelOutput/lmer_Salps_log.rds")
rm(m1_salp)
graphics.off()

########## Jellyfish ########## 
dat_Jelly <-  dat_Jelly %>% filter(Tow != "S3" & Tow != "S2" & Tow != "S1")
# m1_jelly <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
#              ns(Bathy, 3) + log10(Chl) + 
#              exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
#              (1|Gear_Mesh),
#            data = dat_Jelly, weights = WtVec)
m1_jelly <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                   ns(Bathy, 3) + log10(Chl) + 
                   exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                   (1|Gear_Mesh),
                 data = dat_Jelly)
summary(m1_jelly)
anova(m1_jelly) # All Sig
plot(m1_jelly)
fSummariseLMERs(m1_jelly, "Jelly") # r2 = 64.8% (fixed = 12.9%)
fPlotLMERs_Abundance(m1_jelly, "Jellyfish", "log")
saveRDS(m1_jelly, file = "ModelOutput/lmer_Jellyfish_log.rds")
rm(m1_jelly)
graphics.off()
