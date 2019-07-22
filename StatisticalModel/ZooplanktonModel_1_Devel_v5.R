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
# library(arsenal)
library(tidyverse)
library(lmer)
library(visreg)
library(splines)
library(lme4)

source("fHarmonic.R") # Harmonic function
source("fAddWeighting.R")
source("fPlotLMERs.R")

dat <- readRDS("LatestDatabaseOuput_Final_Enviro.rds") # n=985,344

dat <- dat %>% 
  mutate(
    HarmHour = (TimeLocal/24)*2*pi, # Convert to radians
    HarmDOY = (DOY2/365)*2*pi, # Convert to radians
    Latitude2 = abs(Latitude),
    TotAbundance = replace(TotAbundance, TotAbundance > 10000, 10000),
    Bathy = replace(Bathy, Bathy > 6000, 6000),
    Mid_Z = replace(Mid_Z, Mid_Z > 1000, 1000),
    Mesh = replace(Mesh, Mesh > 500, 500), 
    SST = replace(SST, SST > 28, 28)
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
# NOTE: gamm4 and gamm are too slow - will not even run simple models even if increase memory
# NOTE: changing method = to "ML" or "REML" does not fix the random effects

# Test Gear
# dat_Euph$Gear_Msh <- as.factor(paste(as.character(dat_Euph$Gear), as.character(dat_Euph$Mesh), sep =  "_"))
# grepl in fPlotGAMMs picks up Gear_Mesh as Mesh if spell it in full

ggplot(data = dat_Euph, aes(x = log10(Chl), y = log10(TotAbundance + Min))) + 
  geom_point() + geom_text(aes(label = Project))
# BCF - POFI is very high with low Chl-a
dat_Euph <- dat_Euph %>% filter(Project != "BCF - POFI")

# What are the high euph values at high SST?
ggplot(data = dat_Euph, aes(x = SST, y = log10(TotAbundance + Min))) + 
         geom_point() # No visible increase at Equator

ggplot(data = dat_Euph, aes(x = SST, y = DOY2)) + 
  geom_point(aes(colour = log10(TotAbundance + Min))) # Warm temperatures in summer - try cutting off at 30oC

ggplot(data = dat_Euph, aes(x = SST, y = DOY2)) + 
  geom_point(aes(colour = log10(TotAbundance + Min)))

dat_Euph2 <- dat_Euph %>% mutate(SSTbreaks = cut(SST, breaks = -1:28))

ggplot(data = dat_Euph2, aes(x = DOY2, y = log10(TotAbundance + Min))) + 
  facet_wrap(vars(SSTbreaks)) +
  geom_point() + 
  geom_smooth(method = "loess")

summary(dat_Euph$SST)
dat_Euph <- dat_Euph %>% mutate(SST = replace(SST, SST > 28, 28))

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
with(dat, cor(Bathy, Chl))
dat_Euph <- dat_Euph %>% drop_na(Gear)

# Even with lmer, using 1|Gear leads to euphs going up at Equator in winter
m1 <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
           Tow + Mesh + 
           ns(Bathy, 3) + log(Chl) + 
           exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
         data = dat_Omni, weights = WtVec)

m1 <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
           Tow + Mesh + 
           ns(Bathy, 3) + log(Chl) + 
           exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
         data = dat_Carn, weights = WtVec)

m1 <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
           Tow + Mesh + 
           ns(Bathy, 3) + log(Chl) + 
           exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
         data = dat_Chaet, weights = WtVec)

m1 <- lm(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
           Tow + Mesh + 
           ns(Bathy, 3) + log(Chl) + 
           exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
         data = dat_Larv, weights = WtVec)


m1 <- lm(log10(TotAbundance + Min) ~ #te(SST, DOY2, k = 4, bs = c("cr","cc")) +
           fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + ns(Bathy, 3) + Tow + Mesh +
           exp(Chl) + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1),
         data = dat_Euph, weights = WtVec)


m1 <- lm(log10(TotAbundance + Min) ~ #te(SST, DOY2, k = 4, bs = c("cr","cc")) +
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + ns(Bathy, 3) + Mesh + Type +
             exp(Chl) + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2),
           data = dat_Euph, weights = WtVec)

m1 <- lmer(log10(TotAbundance + Min) ~ #te(SST, DOY2, k = 4, bs = c("cr","cc")) +
            fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + ns(Bathy, 3) + Mesh + Tow +
            log10(Chl) + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + (1|Transect),
          data = dat_Omni, weights = WtVec)
summary(m1)
anova(m1)

m1 <- gam(log10(TotAbundance + Min) ~ #te(SST, DOY2, k = 4, bs = c("cr","cc")) +
            fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + s(Bathy, k = 3) + 
            exp(Chl) + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) + 
            s(Gear, bs = "re"),
          data = dat_Euph, weights = WtVec)

m1 <- gam(log10(TotAbundance + Min) ~ te(SST, DOY2, k = 2, bs = c("cr","cc")) +
#                fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + s(Bathy, k = 3) + 
                Mesh + log(Chl) + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
                s(Transect, bs = "re"),
               data = dat_Larv, weights = WtVec)


library(colorRamps) # for Matlab like colour scheme

### Omni ###
m1 <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) + (1|Longhurst),
           data = dat_Omni, weights = WtVec)
x11(width = 15, height = 6, title = "Omni")
par(mfrow = c(2,5), mar = c(4,4,2,2))
# NOTE: To get CIs for lmer, need to use type = contrast (random effects cancel)
visreg2d(m1, yvar = "HarmDOY", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
         ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "deepskyblue2")
visreg(m1, "SST", type = "contrast", scale = "response", band = TRUE, rug = F)
visreg(m1, "HarmDOY", type = "contrast", scale = "response")
visreg(m1, "Chl", type = "contrast", scale = "response")
visreg(m1, "Bathy", type = "contrast", scale = "response")
visreg(m1, "Mesh", type = "contrast", scale = "response")
visreg(m1, "Tow", type = "contrast", scale = "response")
visreg(m1, "Mid_Z", type = "contrast", scale = "response", xlab = "Depth")
visreg(m1, "HarmHour", type = "contrast", scale = "response", xlab = "Time of Day")
visreg(m1, "HarmHour", by = "Mid_Z", type = "contrast", scale = "response", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))
# Bands for HarmHour by Mid_Z come out on pdf not eps, but surface doesn't come out on pdf properly
dev.copy2eps(file = "Omni_LMER.eps")

# Carn
m1 <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow + ns(Bathy, 3) 
             + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) + (1|Longhurst),
           data = dat_Carn, weights = WtVec)
x11(width = 15, height = 6, title = "Carn")
par(mfrow = c(2,5), mar = c(4,4,2,2))
visreg2d(m1, yvar = "HarmDOY", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
         ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "deepskyblue2")
visreg(m1, "SST", type = "contrast", scale = "response", band = TRUE, rug = F)
visreg(m1, "HarmDOY", type = "contrast", scale = "response")
visreg(m1, "Chl", type = "contrast", scale = "response")
visreg(m1, "Bathy", type = "contrast", scale = "response")
visreg(m1, "Mesh", type = "contrast", scale = "response")
visreg(m1, "Tow", type = "contrast", scale = "response")
visreg(m1, "Mid_Z", type = "contrast", scale = "response", xlab = "Depth")
visreg(m1, "HarmHour", type = "contrast", scale = "response", xlab = "Time of Day")
visreg(m1, "HarmHour", by = "Mid_Z", type = "contrast", scale = "response", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))
dev.copy2eps(file = "Carn_LMER.eps")

# Euph
# Need to use Project rather than Transect otherwise negative relationship with Chl
m1 <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + exp(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Project) + (1|Longhurst),
           data = dat_Euph, weights = WtVec)
x11(width = 15, height = 6, title = "Euph")
par(mfrow = c(2,5), mar = c(4,4,2,2))
visreg2d(m1, yvar = "HarmDOY", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
         ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "deepskyblue2")
visreg(m1, "SST", type = "contrast", scale = "response", band = TRUE, rug = F)
visreg(m1, "HarmDOY", type = "contrast", scale = "response")
visreg(m1, "Chl", type = "contrast", scale = "response")
visreg(m1, "Bathy", type = "contrast", scale = "response")
visreg(m1, "Mesh", type = "contrast", scale = "response")
visreg(m1, "Tow", type = "contrast", scale = "response")
visreg(m1, "Mid_Z", type = "contrast", scale = "response", xlab = "Depth")
visreg(m1, "HarmHour", type = "contrast", scale = "response", xlab = "Time of Day")
visreg(m1, "HarmHour", by = "Mid_Z", type = "contrast", scale = "response", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))
dev.copy2eps(file = "Euph_LMER.eps")

# Larv
m1 <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) +  (1|Longhurst),
           data = dat_Larv, weights = WtVec)
x11(width = 15, height = 6, title = "Larv")
par(mfrow = c(2,5), mar = c(4,4,2,2))
visreg2d(m1, yvar = "HarmDOY", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
         ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "deepskyblue2")
visreg(m1, "SST", type = "contrast", scale = "response", band = TRUE, rug = F)
visreg(m1, "HarmDOY", type = "contrast", scale = "response")
visreg(m1, "Chl", type = "contrast", scale = "response")
visreg(m1, "Bathy", type = "contrast", scale = "response")
visreg(m1, "Mesh", type = "contrast", scale = "response")
visreg(m1, "Tow", type = "contrast", scale = "response")
visreg(m1, "Mid_Z", type = "contrast", scale = "response", xlab = "Depth")
visreg(m1, "HarmHour", type = "contrast", scale = "response", xlab = "Time of Day")
visreg(m1, "HarmHour", by = "Mid_Z", type = "contrast", scale = "response", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))
dev.copy2eps(file = "Larv_LMER.eps")

# Chaet
m1 <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) +  (1|Longhurst),
           data = dat_Chaet, weights = WtVec)
x11(width = 15, height = 6, title = "Chaet")
par(mfrow = c(2,5), mar = c(4,4,2,2))
visreg2d(m1, yvar = "HarmDOY", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
         ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "deepskyblue2")
visreg(m1, "SST", type = "contrast", scale = "response", band = TRUE, rug = F)
visreg(m1, "HarmDOY", type = "contrast", scale = "response")
visreg(m1, "Chl", type = "contrast", scale = "response")
visreg(m1, "Bathy", type = "contrast", scale = "response")
visreg(m1, "Mesh", type = "contrast", scale = "response")
visreg(m1, "Tow", type = "contrast", scale = "response")
visreg(m1, "Mid_Z", type = "contrast", scale = "response", xlab = "Depth")
visreg(m1, "HarmHour", type = "contrast", scale = "response", xlab = "Time of Day")
visreg(m1, "HarmHour", by = "Mid_Z", type = "contrast", scale = "response", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))
dev.copy2eps(file = "Chaet_LMER.eps")

# Salp
m1 <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) + (1|Longhurst),
           data = dat_Salp, weights = WtVec)
x11(width = 15, height = 6, title = "Salp")
par(mfrow = c(2,5), mar = c(4,4,2,2))
visreg2d(m1, yvar = "HarmDOY", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
         ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "deepskyblue2")
visreg(m1, "SST", type = "contrast", scale = "response", band = TRUE, rug = F)
visreg(m1, "HarmDOY", type = "contrast", scale = "response")
visreg(m1, "Chl", type = "contrast", scale = "response")
visreg(m1, "Bathy", type = "contrast", scale = "response")
visreg(m1, "Mesh", type = "contrast", scale = "response")
visreg(m1, "Tow", type = "contrast", scale = "response")
visreg(m1, "Mid_Z", type = "contrast", scale = "response", xlab = "Depth")
visreg(m1, "HarmHour", type = "contrast", scale = "response", xlab = "Time of Day")
visreg(m1, "HarmHour", by = "Mid_Z", type = "contrast", scale = "response", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))
dev.copy2eps(file = "Salp_LMER.eps")

# Jelly
m1 <- lmer(log10(TotAbundance + Min) ~ fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             Mesh + Tow +
             ns(Bathy, 3) + log10(Chl) + 
             exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1) + 
             (1|Transect) + (1|Longhurst),
           data = dat_Jelly, weights = WtVec)
x11(width = 15, height = 6, title = "Jelly")
par(mfrow = c(2,5), mar = c(4,4,2,2))
visreg2d(m1, yvar = "HarmDOY", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
         ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "deepskyblue2")
visreg(m1, "SST", type = "contrast", scale = "response", band = TRUE, rug = F)
visreg(m1, "HarmDOY", type = "contrast", scale = "response")
visreg(m1, "Chl", type = "contrast", scale = "response")
visreg(m1, "Bathy", type = "contrast", scale = "response")
visreg(m1, "Mesh", type = "contrast", scale = "response")
visreg(m1, "Tow", type = "contrast", scale = "response")
visreg(m1, "Mid_Z", type = "contrast", scale = "response", xlab = "Depth")
visreg(m1, "HarmHour", type = "contrast", scale = "response", xlab = "Time of Day")
visreg(m1, "HarmHour", by = "Mid_Z", type = "contrast", scale = "response", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))
dev.copy2eps(file = "Jelly_LMER.eps")


x11(width = 15, height = 6)
par(mfrow = c(2,5), mar = c(4,4,2,2))
visreg2d(m1, yvar = "HarmDOY", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
         ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)")
visreg(m1, "SST", scale = "response")
visreg(m1, "HarmDOY", scale = "response")
visreg(m1, "Chl", scale = "response")
visreg(m1, "Bathy", scale = "response")
visreg(m1, "Mesh", scale = "response")
visreg(m1, "Tow", scale = "response")
visreg(m1, "Mid_Z", scale = "response", xlab = "Depth")
visreg(m1, "HarmHour", scale = "response", xlab = "Time of Day")
visreg(m1, "HarmHour", by = "Mid_Z", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", scale = "response", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))
dev.copy2pdf(file = "Omni_NoFixedOrRandomSamplingEffects.pdf", paper = "A4r")
summary(m1)
              
m1_euph <- gam(log10(TotAbundance + Min) ~ te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 # ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) + Mesh +
                 exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Gear, bs = "re") + s(Transect, bs = "re") + s(Longhurst, bs = "re"),
               data = dat_Euph, weights = WtVec)
summary(m1_euph) # r2 = 44.9 % No transect and Gear as fixed
# r2=44.9%, No transect and Gear as random
# r2=40.8%, No transect, No Longhurst, and Gear as random - Makes Chla negative, but Mesh n.s.
# r2=53.4%, Transect included, No Longhurst, and Gear as random
# r2=54.6%, Transect included, Longhurst included, and Gear as random

# r2=42.9% Gear_Mesh fixed, No transect, No Longhurst
# r2=42.9% Gear_Mesh random, No transect, No Longhurst
# r2=54.3% Gear_Mesh random, Transect random, No Longhurst
# r2=55.4% Gear_Mesh random, Transect random, Longhurst random

# r2=54.6% Gear random, Transect random, Longhurst random, Mesh fixed
# r2=54.6% Gear random, Transect random, Longhurst random, Mesh fixed (SST <= 30oC) - goes up less at Equator
# r2=54.6% Gear random, Transect random, Longhurst random, Mesh fixed (SST <= 28oC) - goes up less at Equator

# Tensor instead of ti() # r2 = 54.5%
# Use tensor instead of ti(). Use 28oC.

anova(m1_euph)
visreg2d(m1_euph, yvar = "DOY2", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
         ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)")
fPlotGAMMs(m1_euph, "Euphausiids")
dev.copy2pdf(file = "Euphausiids_GearRandomTransectLonghurst_MeshFixed_SST30.pdf", paper = "A4r")
visreg(m1_euph, "Gear_Mesh", scale = "response")


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
