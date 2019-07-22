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

source("fHarmonic.R") # Harmonic function
source("fAddWeighting.R") # Because CPR is 90% of the data, use weights for individual Groups to give CPR 50% weight and Nets 50% weight
source("fPlotLMERs_Abundance.R") # Plot the Linear Mixed Effects Models
source("fPlotAbundanceLM.R") # Plot the Linear Mixed Effects Models
source("fPlotAbundanceLM_GearMesh.R") # Plot the Linear Mixed Effects Models
source("fSummariseLMERs.R") # Summaries calculating r2 and diagnostic plots

dat <- readRDS("LatestDatabaseOuput_Final_Enviro.rds") # n = 1,009,208

dat <- dat %>% 
  mutate(
    HarmHour = (TimeLocal/24)*2*pi, # Convert to radians
    HarmDOY = (DOY2/365)*2*pi, # Convert to radians
    Latitude2 = abs(Latitude),
    TotAbundance = replace(TotAbundance, TotAbundance > 10000, 10000),
    Bathy = replace(Bathy, Bathy > 6000, 6000),
    Mid_Z = replace(Mid_Z, Mid_Z > 1000, 1000),
    Mesh = replace(Mesh, Mesh > 500, 500), 
    SST = replace(SST, SST > 28, 28),
    Chl = replace(Chl, Chl > 10, 10) # This has already been set to 10, but doing it again here for completeness and traceability
  ) %>% 
  mutate(
    Tow = case_when(Project == "SO-CPR" ~ "S1",
                    Project == "IMOS-CPR" ~ "S2",
                    Project == "SAHFOS" ~ "S3",
                    Tow == "H" ~ "H",
                    Tow == "V" ~ "V",
                    Tow == "O" ~ "O")) %>%
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
# Gear_Mesh
m1_omni <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                Gear_Mesh + 
                ns(Bathy, 3) + log10(Chl) + 
                exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
              data = dat_Omni, weights = WtVec)
summary(m1_omni)
anova(m1_omni) # All significant
plot(m1_omni)
# SummariseLMERs(m1_carn, "Carnivores") # r2 = 73.7% (fixed = 9.65%) OLD r2 = 70.0% (fixed = 9.4%)
fPlotAbundanceLM_GearMesh(m1_omni, "OmnivorousCopepods", "log")

m1_omni <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) + (1|Longhurst),
           data = dat_Omni, weights = WtVec)
summary(m1_omni)
anova(m1_omni) # All significant
plot(m1_omni)
fSummariseLMERs(m1_omni, "Omnivores") # r2 = 68.5% (fixed = 18.0%) OLD r2 = 69.9% (fixed = 16.9%)
fPlotLMERs_Abundance(m1_omni, "OmnivorousCopepods", "log")
saveRDS(m1_omni, file = "ModelOutput/lmer_OmniCopepods_log.rds")

# No random effects. Longhurst and Transect Chl range = 0.65. No Longhurst Chl range = 0.55. No Transect Chl range = 0.7 (and less range of SST) 
m1_omni <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  Mesh + Tow +
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
                data = dat_Omni, weights = WtVec)
fPlotAbundanceLM(m1_omni, "OmnivorousCopepods", "log")

# sqrt(sqrt())
m1_omni <- lmer(sqrt(sqrt(TotAbundance)) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  Mesh + Tow +
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                  (1|Transect) + (1|Longhurst),
                data = dat_Omni, weights = WtVec)
summary(m1_omni)
anova(m1_omni) # All significant
plot(m1_omni)
fSummariseLMERs(m1_omni, "Omnivores") # r2 = 74.4% (fixed = 15.5%)
fPlotLMERs(m1_omni, "OmnivorousCopepods", "sqrtsqrt")
saveRDS(m1_omni, file = "ModelOutput/lmer_OmniCopepods_sqrtsqrt.rds")
rm(m1_omni)
graphics.off()


########## CarnCopepods ##########
m1_carn <- lm(log10(TotAbundance + Min) ~ log10(Chl),
              data = dat_Carn, weights = WtVec) # Chl always -ve
summary(m1_carn)

m1_salp <- lm(log10(TotAbundance + Min) ~ log10(Chl),
              data = dat_Salp, weights = WtVec)
summary(m1_salp)

# Gear_Mesh
m1_carn <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  Gear_Mesh + 
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
                data = dat_Carn, weights = WtVec)
summary(m1_carn)
anova(m1_carn) # All significant
plot(m1_carn)
fSummariseLMERs(m1_carn, "Carnivores") # r2 = 73.7% (fixed = 9.65%) OLD r2 = 70.0% (fixed = 9.4%)
fPlotAbundanceLM_GearMesh(m1_carn, "CarnivorousCopepods", "log")
saveRDS(m1_carn, file = "ModelOutput/lmer_CarnCopepods_log.rds")

m1_carn <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow + 
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) + (1|Longhurst),
           data = dat_Carn, weights = WtVec)
summary(m1_carn)
anova(m1_carn) # All significant
plot(m1_carn)
fSummariseLMERs(m1_carn, "Carnivores") # r2 = 73.7% (fixed = 9.65%) OLD r2 = 70.0% (fixed = 9.4%)
fPlotLMERs(m1_carn, "CarnivorousCopepods", "log")
saveRDS(m1_carn, file = "ModelOutput/lmer_CarnCopepods_log.rds")

# No random effects
m1_carn <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  Mesh + Tow + 
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
                data = dat_Carn, weights = WtVec)
fPlotAbundanceLM(m1_carn, "CarnivorousCopepods", "log")

# Gear - fixed
m1_carn <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                Mesh + Tow + 
                ns(Bathy, 3) + log10(Chl) + 
                exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + Gear,
              data = dat_Carn, weights = WtVec)
fPlotAbundanceLM(m1_carn, "CarnivorousCopepods_Gear", "log")
summary(m1_carn)

# Gear - random
m1_carn <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                Mesh +  
                ns(Bathy, 3) + log10(Chl) + 
                exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + (1|Gear),
              data = dat_Carn, weights = WtVec)
fPlotLMERs_Abundance(m1_carn, "CarnivorousCopepods", "log_GearRandom")
summary(m1_carn)

# Just Chl - strongly negative
m1_carn <- lmer(log10(TotAbundance + Min) ~ log10(Chl) + (1|Gear_Mesh),
                data = dat_Carn, weights = WtVec)
summary(m1_carn)
fPlotAbundanceLM(m1_carn, "CarnivorousCopepods_ChlOnly", "log")

m1_carn <- lmer(sqrt(sqrt(TotAbundance)) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  Mesh + Tow + 
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                  (1|Transect) + (1|Longhurst),
                data = dat_Carn, weights = WtVec)
summary(m1_carn)
anova(m1_carn) # All significant
plot(m1_carn)
fSummariseLMERs(m1_carn, "Carnivores") # r2 = 74.1% (fixed = 16.7%)
fPlotLMERs(m1_carn, "CarnivorousCopepods", "sqrtsqrt")
saveRDS(m1_carn, file = "ModelOutput/lmer_CarnCopepods_sqrtsqrt.rds")
rm(m1_carn)
graphics.off()


########## Chaetognaths ##########
m1_chaet <- lm(log10(TotAbundance + Min) ~ log10(Chl),
              data = dat_Chaet, weights = WtVec)
summary(m1_chaet)

# Gear_Mesh
m1_chaet <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                Gear_Mesh + 
                ns(Bathy, 3) + log10(Chl) + 
                exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
              data = dat_Chaet, weights = WtVec)
summary(m1_chaet)
anova(m1_chaet) # All significant
plot(m1_chaet)
# SummariseLMERs(m1_carn, "Carnivores") # r2 = 73.7% (fixed = 9.65%) OLD r2 = 70.0% (fixed = 9.4%)
fPlotAbundanceLM_GearMesh(m1_chaet, "Chaetognaths", "log")

m1_chaet <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) +  (1|Longhurst),
           data = dat_Chaet, weights = WtVec)
summary(m1_chaet)
anova(m1_chaet) # All significant
plot(m1_chaet)
fSummariseLMERs(m1_chaet, "Chaetognaths") #r2 = 61.0% (fixed = 6.1%) OLD r2 = 62.0% (fixed = 6.6%) 
fPlotLMERs(m1_chaet, "Chaetognaths", "log") 
saveRDS(m1_chaet, file = "ModelOutput/lmer_Chaetognaths_log.rds")

# No random effects
m1_chaet <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                   Mesh + Tow +
                   ns(Bathy, 3) + log10(Chl) + 
                   exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
                 data = dat_Chaet, weights = WtVec)
fPlotAbundanceLM(m1_chaet, "Chaetognaths", "log")
plot(m1_chaet)

# Gear
m1_chaets <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                Mesh + Tow +
                ns(Bathy, 3) + log10(Chl) + 
                exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + Gear,
              data = dat_Chaet, weights = WtVec)
summary(m1_chaets)
fPlotAbundanceLM(m1_chaets, "Chaetognaths_Gear", "log")

m1_chaet <- lmer(sqrt(sqrt(TotAbundance)) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                   Mesh + Tow +
                   ns(Bathy, 3) + log10(Chl) + 
                   exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                   (1|Transect) +  (1|Longhurst),
                 data = dat_Chaet, weights = WtVec)
summary(m1_chaet)
anova(m1_chaet) # All significant
plot(m1_chaet)
fSummariseLMERs(m1_chaet, "Chaetognaths") # r2 = 61.7% (fixed = 14.5%)
fPlotLMERs(m1_chaet, "Chaetognaths", "sqrtsqrt")
saveRDS(m1_chaet, file = "ModelOutput/lmer_Chaetognaths_sqrtsqrt.rds")
rm(m1_chaet)
graphics.off()


########## Larvaceans ##########
m1_larv <- lm(log10(TotAbundance + Min) ~ log10(Chl),
              data = dat_Larv, weights = WtVec)
summary(m1_larv)

# Gear_Mesh
m1_larv <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                Gear_Mesh + 
                ns(Bathy, 3) + log10(Chl) + 
                exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
              data = dat_Larv, weights = WtVec)
summary(m1_larv)
anova(m1_larv) # All significant
plot(m1_larv)
# SummariseLMERs(m1_carn, "Carnivores") # r2 = 73.7% (fixed = 9.65%) OLD r2 = 70.0% (fixed = 9.4%)
fPlotAbundanceLM_GearMesh(m1_larv, "Larvaceans", "log")

m1_larv <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) +  (1|Longhurst),
           data = dat_Larv, weights = WtVec)
summary(m1_larv)
anova(m1_larv) # n.s. exp(-Mid_Z/1000):fHarmonic(HarmHour, k = 1)

m1_larv <- update(m1_larv, . ~ . -exp(-Mid_Z/1000):fHarmonic(HarmHour, k = 1))
anova(m1_larv) # n.s. fHarmonic(HarmHour, k = 1)

m1_larv <- update(m1_larv, . ~ . -fHarmonic(HarmHour, k = 1))
anova(m1_larv) # All sig. now
plot(m1_larv)

fSummariseLMERs(m1_larv, "Larvaceans") # r2 = 68.8% (fixed = 15.6%) OLD r2 = 70.0% (fixed = 9.4%)
fPlotLMERs(m1_larv, "Larvaceans", "log")
saveRDS(m1_larv, file = "ModelOutput/lmer_Larvaceans_log.rds")

# No random effects
m1_larv <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  Mesh + Tow +
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
                data = dat_Larv, weights = WtVec)
fPlotAbundanceLM(m1_larv, "Larvaceans", "log")
plot(m1_larv)

# Gear
m1_larv <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                Mesh + Tow +
                ns(Bathy, 3) + log10(Chl) + 
                exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + Gear,
              data = dat_Larv, weights = WtVec)
summary(m1_larv)
fPlotAbundanceLM(m1_larv, "Larvaceans_Gear", "log")



m1_larv <- lmer(sqrt(sqrt(TotAbundance)) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  Mesh + Tow + 
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                  (1|Transect) + (1|Longhurst),
                data = dat_Carn, weights = WtVec)
summary(m1_larv)
anova(m1_larv) # All sig.
plot(m1_larv)
fSummariseLMERs(m1_larv, "Larvaceans") # r2 = 74.1% (fixed = 16.7%)
fPlotLMERs(m1_larv, "Larvaceans", "sqrtsqrt")
saveRDS(m1_larv, file = "ModelOutput/lmer_Larvaceans_sqrtsqrt.rds")
rm(m1_larv)
graphics.off()


########## Euphausiids ##########
# CPR data dodgy - says 97% no euphs - wrong by Richardson et al. (2006) - should be about 60%
# NOTE: gamm4 and gamm are too slow - will not even run simple models even if increase memory
# NOTE: changing method = to "ML" or "REML" does not fix the random effects
# Type is good, but flattens it out a bit
# Gear pulls it up at Equator
# Gear_Mesh pulls it up at Equator
# s(Gear_Mesh) pulls it up at Equator
# Mesh flattens it out at Equator
# s(Transect) changes it at the Pole and at Equator (pulls up)
# s(Project) changes it even more at the Pole and at Equator (pulls up)
# s(Longhurst) changes the shape a fair bit too, but better than s(Transect) and s(Project)
# These are same as fixed effects
# Type + Gear - kicks up
# Tow changes it a lot!
# Bathy can lessen effect of Chl-a (and Chl-a is strongly positive with all just with surface - except weakly -ve for Chaets)
# Chaets positive with Chl-a if SST is ignored...
# Chaets positive with Chl-a if SST not in as surface - just as s(SST)...
# For euphausiids, log10(Chl) is n.s., but exp(Chl) is massively significant...
# Choose from linear, exp or log10 for each species

# Gear_Mesh
m1_euph <- lm(log10(TotAbundance + Min) ~ log10(Chl),
              data = dat_Euph, weights = WtVec)
summary(m1_euph)

m1_euph <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                Gear_Mesh + 
                ns(Bathy, 3) + exp(Chl) + 
                exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
              data = dat_Euph, weights = WtVec)
summary(m1_euph)
anova(m1_euph) # All significant
plot(m1_euph)
# SummariseLMERs(m1_carn, "Carnivores") # r2 = 73.7% (fixed = 9.65%) OLD r2 = 70.0% (fixed = 9.4%)
fPlotAbundanceLM_GearMesh(m1_euph, "Euphausiids", "log")

# Need to use Project rather than Transect otherwise negative relationship with Chl
m1_euph <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + exp(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Project) + (1|Longhurst),
           data = dat_Euph, weights = WtVec)
summary(m1_euph)
anova(m1_euph) # All sig.
plot(m1_euph)
fSummariseLMERs(m1_euph, "Euphausiids") # r2 = 72.5% (fixed = 4.5%) OLD: r2 = 72.5% (fixed = 5.5%)
fPlotLMERs(m1_euph, "Euphausiids", "log")
saveRDS(m1_euph, file = "ModelOutput/lmer_Euphausiids_log.rds")

m1_euph <- lmer(sqrt(sqrt(TotAbundance)) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  Mesh + Tow +
                  ns(Bathy, 3) + exp(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                  (1|Project) + (1|Longhurst),
                data = dat_Euph, weights = WtVec)
summary(m1_euph)
anova(m1_euph) # All sig.
plot(m1_euph)
fSummariseLMERs(m1_euph, "Euphausiids") # r2 = 77.7% (fixed = 10.4%)
fPlotLMERs(m1_euph, "Euphausiids", "sqrtsqrt")
saveRDS(m1_euph, file = "ModelOutput/lmer_Euphausiids_sqrt.rds")
rm(m1_euph)
graphics.off()



########## Salps ########## 
m1_salp <- lm(log10(TotAbundance + Min) ~ log10(Chl),
              data = dat_Salp, weights = WtVec)
summary(m1_salp)

# Gear_Mesh
m1_salp <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                Gear_Mesh + 
                ns(Bathy, 3) + log10(Chl) + 
                exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
              data = dat_Salp, weights = WtVec)
summary(m1_salp)
anova(m1_salp) # All significant
plot(m1_salp)
# SummariseLMERs(m1_carn, "Carnivores") # r2 = 73.7% (fixed = 9.65%) OLD r2 = 70.0% (fixed = 9.4%)
fPlotAbundanceLM_GearMesh(m1_salp, "Euphausiids", "log")

m1_salp <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) + (1|Longhurst),
           data = dat_Salp, weights = WtVec)
summary(m1_salp)
anova(m1_salp) # Mesh now p =0.09, but still +ve OLD: Mesh n.s. (p = 0.06), but Mesh positive

m1_salp <- update(m1_salp, . ~ . -Mesh)
anova(m1_salp) # All sig. now
plot(m1_salp)
fSummariseLMERs(m1_salp, "Salps") # r2 = 74.7% (fixed = 12.4%) OLD: r2 = 70.9% (fixed = 13.2%)
fPlotLMERs(m1_salp, "Salps", "log")
saveRDS(m1_salp, file = "ModelOutput/lmer_Salps_log.rds")

m1_salp <- lmer(sqrt(sqrt(TotAbundance)) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                  Mesh + Tow +
                  ns(Bathy, 3) + log10(Chl) + 
                  exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                  (1|Transect) + (1|Longhurst),
                data = dat_Salp, weights = WtVec)
summary(m1_salp)
anova(m1_salp) # Mesh n.s. (p = 0.35), but Mesh positive

m1_salp <- update(m1_salp, . ~ . -Mesh)
anova(m1_salp) # All sig. now
plot(m1_salp)
fSummariseLMERs(m1_salp, "Salps") # r2 = 56.8% (fixed = 4.5%)
fPlotLMERs(m1_salp, "Salps", "sqrtsqrt")
saveRDS(m1_salp, file = "ModelOutput/lmer_Salps_sqrtsqrt.rds")
rm(m1_salp)
graphics.off()



########## Jellyfish ########## 
m1_jelly <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) + (1|Longhurst),
           data = dat_Jelly, weights = WtVec)
summary(m1_jelly)
anova(m1_jelly) # All Sig
plot(m1_jelly)
fSummariseLMERs(m1_jelly, "Jelly") # r2 = 71.6% (fixed = 03.2%) OLD: r2 = 62.3% (fixed = 13.7%)
fPlotLMERs(m1_jelly, "Jellyfish", "log")
saveRDS(m1_jelly, file = "ModelOutput/lmer_Jellyfish_log.rds")

m1_jelly <- lmer(sqrt(sqrt(TotAbundance)) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
                   Mesh + Tow +
                   ns(Bathy, 3) + log10(Chl) + 
                   exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                   (1|Transect) + (1|Longhurst),
                 data = dat_Jelly, weights = WtVec)
summary(m1_jelly)
anova(m1_jelly) # n.s. exp(-Mid_Z/1000):fHarmonic(HarmHour, k = 1)
m1_jelly <- update(m1_jelly, . ~ . -exp(-Mid_Z/1000):fHarmonic(HarmHour, k = 1))
anova(m1_jelly) # n.s. All sig. now
plot(m1_jelly)
fSummariseLMERs(m1_jelly, "Jelly") # r2 = 62.3% (fixed = 13.7%)
fPlotLMERs(m1_jelly, "Jellyfish", "sqrtsqrt")
saveRDS(m1_jelly, file = "ModelOutput/lmer_Jellyfish_sqrtsqrt.rds")
rm(m1_jelly)
graphics.off()

