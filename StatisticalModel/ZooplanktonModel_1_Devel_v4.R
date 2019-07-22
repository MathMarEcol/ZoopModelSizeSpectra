# ZooplanktonModelDevel.R
# Ant, Jase, Ryan
# Created 11th April 2019
# Last Edited: 14th June 2019

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
# https://stats.stackexchange.com/questions/197952/two-methods-of-adding-random-effects-to-a-gam-give-very-different-results-why-i
# Consider a Copepod only model (joining, carn and omni)

############################  Preliminaries ############################  
library(arsenal)
# library(effects)
library(tidyverse)
#@ library(lme4)
# library(MuMIn)
# library(splines)
library(mgcv)
library(visreg)

source("fHarmonic.R") # Harmonic function
source("fAddWeighting.R")
source("fAddWeighting2.R")
source("fPlotGAMMs.R")

dat <- readRDS("LatestDatabaseOuput_Final_Enviro.rds") # n=985,344

dat <- dat %>% 
  mutate(
    HarmHour = (TimeLocal/24)*2*pi, # Convert to radians
    HarmDOY = (DOY2/365)*2*pi, # Convert to radians
    Latitude2 = abs(Latitude),
    TotAbundance = replace(TotAbundance, TotAbundance > 10000, 10000),
    Bathy = replace(Bathy, Bathy > 6000, 6000),
    Mid_Z = replace(Mid_Z, Mid_Z > 1000, 1000),
    Mesh = replace(Mesh, Mesh > 500, 500)
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

dat_Omni <- fAddWeighting(dat, "OmniCopepods")
dat_Carn <- fAddWeighting(dat, "CarnCopepods")
dat_Chaet <- fAddWeighting(dat, "Chaetognaths")
dat_Euph <- fAddWeighting(dat, "Euphausiids")
dat_Jelly <- fAddWeighting(dat, "Jellyfish")
dat_Larv <- fAddWeighting(dat, "Larvaceans")
dat_Salp <- fAddWeighting(dat, "Salps")

########## Omnivorous Copepods ##########
# For the tensor surface the different basis functions ("cr" = cubic regression, 
# "cs" = cubic regression with shrinkage, "tp" = thin plate spline, "ts" = tprs with shrinkage)
# all give the same results. Use "tp", as a bit better when gaps in data
# For s(Bathy), "tp" and "ts" have tighter SEs than "cr" abnd "cs"
# Bathy seems best with s(Bathy, k = 4)
# ti(SST, DOY2, k = 4, bs = c("tp","cc")) - k=4 appears best (less wiggly), but check
# Including Latitude with SST produces Latitude plot that doesn't really make sense (and gives very large SEs)
# s(SST, DOY2) gives a very folded surface - no good
# fHarmonic(k=2) seems generally to be best
# NOTE: fHarmonic(DOY2 * pi / 180, k = 2) * fHarmonic(HarmHour, k = 2) is not useful
# Set up weighting vectors
# Each group has different weighting vector because each has different number of Net and CPR samples

m1_omni <- gam(log10(TotAbundance + Min) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"), 
               weights = WtVec, data = dat_Omni)
# Time = 1 min 28 secs
m1_omni2 <- bam(log10(TotAbundance + Min) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"), 
               weights = WtVec, data = dat_Omni, chunk.size = 75000)
# Time = 41 secs

summary(m1_omni) #r2=67.7%
anova(m1_omni) # All significant
graphics.off()
par(mfrow = c(2,2))
gam.check(m1_omni) # Looks reasonable
fPlotGAMMs(m1_omni2, "OmniCopepods")
dev.copy2pdf(file = "OmniCopepods.pdf", paper = "A4r") #

saveRDS(m1_omni, file = "ModelOutput/gamm_OmniCopepods.rds")
rm(m1_omni)

########## CarnCopepods ##########
m1_carn <- gam(log10(TotAbundance + Min) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"), 
               weights = WtVec, data = dat_Carn)
summary(m1_carn) # r2 = 59.6%
anova(m1_carn) # All sig
fPlotGAMMs(m1_carn, "CarnCopepods")
dev.copy2pdf(file = "CarnCopepods.pdf", paper = "A4r")
saveRDS(m1_carn, file = "ModelOutput/gamm_CarnCopepods.rds")
rm(m1_carn)

########################################
########## Chaetognaths ##########
m1_chaet <- gam(log10(TotAbundance + Min) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) +
                  Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                  s(Longhurst, bs = "re") + s(Transect, bs = "re"),
                data = dat_Chaet, weights = WtVec)
summary(m1_chaet) # r2 = 57.5%
anova(m1_chaet) # All sig
fPlotGAMMs(m1_chaet, "Chaetognaths")
dev.copy2pdf(file = "Chaetognaths.pdf", paper = "A4r")

saveRDS(m1_chaet, file = "ModelOutput/gamm_Chaetognaths.rds")
rm(m1_chaet)

########## Larvaceans ##########
m1_larv <- gam(log10(TotAbundance + Min) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"),
               data = dat_Larv, weights = WtVec)
summary(m1_larv) # r2 = 59.6%
anova(m1_larv) # Interaction exp(-Mid_Z/1000):fHarmonic(HarmHour, k = 2) n.s., so drop

m1_larv <- update(m1_larv, . ~ . -exp(-Mid_Z/1000):fHarmonic(HarmHour, k = 2))
summary(m1_larv) # r2 = 59.6%
anova(m1_larv) # fHarmonic(HarmHour, k = 2) p = 0.0965, so drop

m1_larv <- update(m1_larv, . ~ . -fHarmonic(HarmHour, k = 2))
summary(m1_larv) # r2 = 59.6%
anova(m1_larv) # All sig

fPlotGAMMs(m1_larv, "Larvaceans") 
dev.copy2pdf(file = "Larvaceans.pdf", paper = "A4r")

saveRDS(m1_larv, file = "ModelOutput/gamm_Larvaceans.rds")
write2pdf(summary(m1_larv), "LarvaceanSummary.pdf", quiet = FALSE, title = "Larvaceans") # This will only work if you have tex installed
rm(m1_larv)

########## Euphausiids ##########
# CPR data dodgy - says 97% no euphs - wrong by Richardson et al. (2006) - should be about 60%
m1_euph <- gam(log10(TotAbundance + Min) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"),
               data = dat_Euph, weights = WtVec)
summary(m1_euph) # r2 = 54.9 %
anova(m1_euph) # All sig

graphics.off()
par(mfrow = c(2,2))
gam.check(m1_euph) # Looks reasonable

fPlotGAMMs(m1_euph, "Euphausiids")
dev.copy2pdf(file = "Euphausiids.pdf", paper = "A4r")

saveRDS(m1_euph, file = "ModelOutput/gamm_Euphausiids.rds")
rm(m1_euph)


########## Salps ########## 
ggplot(data = dat_Salp, aes(Mesh, log10(TotAbundance + Min))) +
  geom_point() + geom_smooth(method = "lm")

ggplot(data = dat_Salp, aes(Mesh, log10(TotAbundance + Min))) +
  geom_point() + geom_smooth(method = "lm") + facet_wrap(vars(Tow))
# still +ve slope for Mesh when >350 um max

# Presence/Absence GAM
dat_Salp2 <- dat_Salp %>% mutate(TotAbundance = replace(TotAbundance, TotAbundance > 0, 1))

m1_salps2 <- bam(TotAbundance ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) +
                  Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                  s(Longhurst, bs = "re") + s(Transect, bs = "re"),
                data = dat_Salp2, family = binomial, chunk.size = 150000, weights = WtVec)
summary(m1_salps2) # r2 = 49.9%
anova(m1_salps2)
fPlotGAMMs(m1_salps2, "Salps")
dev.copy2pdf(file = "Salps_BAM.pdf", paper = "A4r")
# BAM better for salps than GAM

m1_salps <- gam(log10(TotAbundance + Min) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) +
                  Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                  s(Longhurst, bs = "re") + s(Transect, bs = "re"),
                data = dat_Salp, weights = WtVec)
summary(m1_salps) # r2 = 62.4 %
anova(m1_salps) # All sig
fPlotGAMMs(m1_salps, "Salps")
dev.copy2pdf(file = "Salps.pdf", paper = "A4r")

saveRDS(m1_salps, file = "ModelOutput/gamm_Salps.rds")
rm(m1_salps, m1_salps2)

########## Jellyfish ########## 
ggplot(data = dat_Jelly, aes(Mesh, log10(TotAbundance + Min))) +
  geom_point() + geom_smooth(method = "lm")
ggplot(data = dat_Jelly, aes(Mesh, log10(TotAbundance + Min))) +
  geom_point() + geom_smooth(method = "lm") + facet_wrap(vars(Tow))

# Presence/Absence GAM
dat_Jelly2 <- dat_Jelly %>% mutate(TotAbundance = replace(TotAbundance, TotAbundance > 0, 1))

m1_jelly <- bam(TotAbundance ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) + 
                  Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                  s(Longhurst, bs = "re") + s(Transect, bs = "re"), family = binomial, chunk.size = 75000,
                data = dat_Jelly2)
# Takes 9 mins to run, but gam doesn't fit in 30 mins!
# GAM better for Jellies than BAM
summary(m1_jelly)
anova(m1_jelly)
# Simon Wood:  If you must have p-values, then anova is better for any model containing factor variables

m1_jelly <- gam(log10(TotAbundance + Min) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) + 
                  Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                  s(Longhurst, bs = "re") + s(Transect, bs = "re"),
                data = dat_Jelly, weights = WtVec)
summary(m1_jelly) # r2 = 61.5 %
anova(m1_jelly) # All Sig

fPlotGAMMs(m1_jelly, "Jellyfish")
dev.copy2pdf(file = "Jellyfish_BAM.pdf", paper = "A4r")

saveRDS(m1_jelly, file = "ModelOutput/gamm_Jellyfish.rds")
rm(m1_jelly, m1_jelly2)
