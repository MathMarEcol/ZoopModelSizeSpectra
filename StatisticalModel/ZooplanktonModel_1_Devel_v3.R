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
library(effects)
library(tidyverse)
library(lme4)
library(MuMIn)
library(splines)
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
  # mutate(
  #   Tow = case_when(Project == "SO-CPR" ~ "S1",
  #                   Project == "IMOS-CPR" ~ "S2",
  #                   Project == "SAHFOS" ~ "S3",
  #                   Tow == "H" ~ "H",
  #                   Tow == "V" ~ "V",
  #                   Tow == "O" ~ "O")
  #   ) %>%
  mutate(
    Tow = case_when(Project == "SO-CPR" ~ "S",
                     Project == "IMOS-CPR" ~ "S",
                     Project == "SAHFOS" ~ "S",
                     Tow == "H" ~ "H",
                     Tow == "V" ~ "V",
                     Tow == "O" ~ "O")
  ) %>%
  droplevels()

############################  GAMs ############################  
# RANDOM EFFECTS - Use best model here and add Longhurst and Gear. Maybe ShipCruise
# Gear as a random effect would be better as we don't need to choose a level in the maps
# HOW does the best model here (G1), work with the other functional groups

# Define min_val for GAM for each group - not worth doing it more elegantly!
Min <- dat %>% 
  filter(TotAbundance > 0) %>% 
  # group_by(Group) %>% 
  group_by(Group, Tow) %>% 
  # summarise(Min = min(TotAbundance) / 2)
  summarise(Min = min(TotAbundance))
Min <- Min %>% unite(col = MinGrp, Group, Tow)

MinChaet <- as.numeric(Min[1, 2])
MinCarn <- as.numeric(Min[2, 2])
MinOmni <- as.numeric(Min[3, 2])
MinEuph <- as.numeric(Min[4, 2])

MinJelly <- as.numeric(Min[5, 2])
MinLarv <- as.numeric(Min[6, 2])
MinSalp <- as.numeric(Min[7,2])

dat_Omni <- fAddWeighting(dat, "OmniCopepods")
dat_Carn <- fAddWeighting(dat, "CarnCopepods")
dat_Chaet <- fAddWeighting(dat, "Chaetognaths")
dat_Euph <- fAddWeighting(dat, "Euphausiids")
dat_Jelly <- fAddWeighting(dat, "Jellyfish")
dat_Larv <- fAddWeighting(dat, "Larvaceans")
dat_Salp <- fAddWeighting(dat, "Salps")

# capture_tally <- dat_Jelly %>% group_by(Tow, Transect) %>% 
#   summarise(Count_of_Captures = n())  
# uni <- unique(capture_tally$Transect)
# 
# 1-length(uni)/dim(capture_tally)[1]
# 
# capture_tally <- dat_Omni %>% group_by(Tow, Transect) %>% 
#   summarise(Count_of_Captures = n())  
# uni <- unique(capture_tally$Transect)
# 1-length(uni)/dim(capture_tally)[1]

# ### Look at data distribution
# ggplot(data = dat, aes(x = log10(TotAbundance + 1))) +
#   geom_histogram() +
#   theme_bw()
# summary(dat$TotAbundance) # No MD
# 
# ggplot(data = dat %>% filter(Type == "CPR"), aes(x = log10(TotAbundance+1))) +
#   geom_histogram() +
#   theme_bw()
# 
# ggplot(data = dat %>% filter(Type == "Net"), aes(x = log10(TotAbundance+1))) +
#   geom_histogram() +
#   theme_bw()

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

m1_omni <- gam(log10(TotAbundance + MinOmni) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) + 
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"), 
               weights = WtVec, data = dat_Omni)

m1_omni <- gam(log10(TotAbundance + MinOmni) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) + 
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Project, bs = "re"), 
               weights = WtVec, data = dat_Omni)
m1_omni <- gam(sqrt(sqrt(TotAbundance)) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 # ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) + 
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"), 
               weights = WtVec, data = dat_Omni)

m1_omni <- gam(log10(TotAbundance + MinOmni) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 # ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 # ns(SST, df = 3) * fHarmonic(HarmDOY, k = 2) +
                 ns(Latitude2, df = 5) * fHarmonic(HarmDOY, k = 2) +
                 ns(SST, df = 3) +
                 fHarmonic(HarmDOY, k = 2) * fHarmonic(HarmHour, k = 2) + 
                 log10(Chl) + s(Bathy, k = 3) + 
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"), 
               weights = WtVec, data = dat_Omni)
visreg2d(m1_omni, yvar = "Latitude2", xvar = "HarmDOY", scale = "response", ylim = c(0, 60))
visreg2d(m1_omni, xvar = "Latitude2", yvar = "HarmDOY", scale = "response", plot.type="persp")
visreg2d(m1_omni, xvar = "HarmHour", yvar = "HarmDOY", scale = "response", plot.type="persp")
visreg2d(m1_omni, yvar = "HarmHour", xvar = "HarmDOY", scale = "response")
visreg2d(m1_omni, yvar = "HarmHour", xvar = "HarmDOY", scale = "response")
visreg(m1_omni, "HarmDOY", by = "Latitude2", breaks = c(0, 40, 60), xlab = "Time of day", ylab = "log10(Abundance)", 
       scale = "response", overlay = TRUE, rug = 0, strip.names = c("0","30", "60"))




# Surface again, but including SST
m1_omni <- gam(log10(TotAbundance + MinOmni) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(Latitude2, bs = "cr") + ti(DOY2, bs = "cc") + ti(Latitude2, DOY2, k = 4, bs = c("cr","cc")) +
                 # ns(SST, df = 3) * fHarmonic(HarmDOY, k = 2) +
                 ns(SST, df = 3) +
                 log10(Chl) + s(Bathy, k = 3) + 
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"), 
               weights = WtVec, data = dat_Omni)


summary(m1_omni) #r2 = 68.3%, log10(Y+MinOmni) w Transect; log with Project r2=65.5%; r2 = 73.7%, sqrt(Y); r2 = 73.5%, sqrt(sqrt(Y)) with Project r2 = 70.9%
anova(m1_omni) # All significant
gam.check(m1_omni) # Looks reasonable I think
fPlotGAMMs(m1_omni, "OmniCopepods")
dev.copy2pdf(file = "OmniCopepods_sqrtsqrt_Project.pdf", paper = "A4r") #

saveRDS(m1_omni, file = "ModelOutput/gamm_OmniCopepods.rds")
rm(m1_omni)

# visreg(m1_omni, "HarmDOY", by = "SST", breaks = c(5, 15, 25), xlab = "Day of Year", ylab = "log10(Abundance)", 
#         scale = "response", overlay = TRUE, rug = 0, strip.names = c("SST=5","SST=15", "SST=25"))
x11(width = 16, height = 6)
par(mfrow = c(2,6), mar = c(4,4,2,2))
visreg(m1_omni, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(m1_omni, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(m1_omni, "Mesh", scale = "response", xlab = "Mesh  (microns)")
visreg(m1_omni, "Tow", rug = FALSE, scale = "response", xlab = "Tow")
visreg(m1_omni, "SST", scale = "response")
visreg2d(m1_omni, yvar = "DOY2", xvar = "Latitude2", scale = "response", plot.type = "persp", theta = -235, phi = -10, r = 100,
         ticktype = "detailed", xlab = "\nLatitude (º)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)")
visreg2d(m1_omni, xvar = "DOY2", yvar = "Latitude2", scale = "response")
visreg(m1_omni, "Latitude2", scale = "response", xlab = "Latitude")
visreg(m1_omni, "HarmDOY", scale = "response", xlab = "Day of year")
visreg(m1_omni, "HarmHour", scale = "response", xlab = "Time of Day")
visreg(m1_omni, "Mid_Z", scale = "response", xlab = "Depth")
visreg(m1_omni, "HarmHour", by = "Mid_Z", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", 
       scale = "response", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))

  dev.copy2pdf(file = "OmniCopepods_Test.pdf", paper = "A4r") #

########## CarnCopepods ##########
m1_carn <- gam(log10(TotAbundance + MinCarn) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"), 
               weights = WtVec, data = dat_Carn)
m1_carn <- gam(sqrt(sqrt(TotAbundance)) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"), 
               weights = WtVec, data = dat_Carn)

m1_carn <- gam(log10(TotAbundance + MinCarn) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 # ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 # ns(SST, df = 3) * fHarmonic(HarmDOY, k = 2) +
                 ns(Latitude2, df = 5) * fHarmonic(HarmDOY, k = 2) +
                 ns(SST, df = 3) +
                 fHarmonic(HarmDOY, k = 2) * fHarmonic(HarmHour, k = 2) + 
                 log10(Chl) + s(Bathy, k = 3) + 
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"), 
               weights = WtVec, data = dat_Carn)

summary(m1_carn) # r2 = 66.0%, log10(Y+MinCarn); r2 = 67.2%, sqrt; r2 = 70.3%; sqrt(sqrt())
anova(m1_carn) # All sig
fPlotGAMMs(m1_carn, "CarnCopepods")

visreg(m1_carn, "HarmDOY", by = "Latitude2", breaks = c(10, 40, 60), xlab = "Time of day", ylab = "log10(Abundance)", 
       scale = "response", overlay = TRUE, rug = 0, strip.names = c("10","40", "60"))


x11(width = 12, height = 9)
par(mfrow = c(3,4), mar = c(4,4,2,2))
visreg(m1_carn, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(m1_carn, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(m1_carn, "Mesh", scale = "response", xlab = "Mesh  (microns)")
visreg(m1_carn, "Tow", rug = FALSE, scale = "response", xlab = "Tow")
visreg(m1_carn, "SST", scale = "response")
visreg(m1_carn, "Latitude2", scale = "response", xlab = "Latitude")
visreg(m1_carn, "HarmDOY", scale = "response", xlab = "Day of year")
visreg(m1_carn, "HarmDOY", by = "Latitude2", breaks = c(10, 30, 50, 70), xlab = "Day of Year", ylab = "log10(Abundance)", 
       scale = "response", overlay = TRUE, rug = 0, strip.names = c("Lat=10","Lat=30", "Lat=50", "Lat=70"))
visreg(m1_carn, "HarmHour", scale = "response", xlab = "Time of Day")
visreg(m1_carn, "HarmHour", by = "HarmDOY", breaks = c(81/365*2*pi, 181/365*2*pi, 264/365*2*pi, 355/365*2*pi), xlab = "Time of Day", ylab = "log10(Abundance)", 
       scale = "response", overlay = TRUE, rug = 0, strip.names = c("Spr Equ","Sum Sol", "Aut Equ", "Win Sol"))
visreg(m1_carn, "Mid_Z", scale = "response", xlab = "Depth")
visreg(m1_carn, "HarmHour", by = "Mid_Z", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", 
       scale = "response", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))

visreg2d(m1_carn, xvar = "HarmDOY", yvar = "Latitude2", scale = "response", plot.type = "persp")
visreg2d(m1_carn, xvar = "HarmDOY", yvar = "Latitude2", scale = "response")

dev.copy2pdf(file = "CarnCopepods_Latitude.pdf", paper = "A4r")

saveRDS(m1_carn, file = "ModelOutput/gamm_CarnCopepods.rds")
rm(m1_carn)
# write2pdf(summary(gmm1_carn), "CarnSummary.pdf",quiet = FALSE, title = "Carnivorous Copeods") # This will only work if you have tex installed

########################################
########## Chaetognaths ##########

# One chaetognath mesh size with 3000 um - make it next maximum mesh size of 710 for Chaets
# dat_Chaet <- dat_Chaet %>% 
#   mutate(Mesh = replace(Mesh, Mesh > 710, 710))

m1_chaet <- gam(log10(TotAbundance + MinChaet) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) +
                  Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                  s(Longhurst, bs = "re") + s(Transect, bs = "re"),
                data = dat_Chaet, weights = WtVec)
m1_chaet <- gam(sqrt(sqrt(TotAbundance)) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) +
                  Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                  s(Longhurst, bs = "re") + s(Transect, bs = "re"),
                data = dat_Chaet, weights = WtVec)

visreg(m1_chaet, "DOY2", by = "SST", breaks = c(10, 20, 30), xlab = "Time of day", ylab = "log10(Abundance)", 
       scale = "response", overlay = TRUE, rug = 0, strip.names = c("10","20", "30"))

summary(m1_chaet) # r2 = 57.5%, log; r2 = 52.4%, sqrt; r2 = 56.3% sqrt(sqrt())
anova(m1_chaet) # All sig
fPlotGAMMs(m1_chaet, "Chaetognaths")
dev.copy2pdf(file = "Chaetognaths_sqrtsqrt.pdf", paper = "A4r")

saveRDS(m1_chaet, file = "ModelOutput/gamm_Chaetognaths.rds")
rm(m1_chaet)
# write2pdf(summary(gmm1_chaet), "ChaetSummary.pdf",quiet = FALSE, title = "Chaetognaths") # This will only work if you have tex installed

########## Larvaceans ##########
m1_larv <- gam(log10(TotAbundance + MinLarv) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"),
               data = dat_Larv, weights = WtVec)
m1_larv <- gam(sqrt(sqrt(TotAbundance)) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"),
               data = dat_Larv, weights = WtVec)

summary(m1_larv) # r2 = 37.7%
anova(m1_larv) # Chl n.s.

m1_larv <- update(m1_larv, . ~ . -log10(Chl))
summary(m1_larv) # r2 = 37.7%, log; r2 = 39.4% sqrt(sqrt())
anova(m1_larv) # Larvacean model has more significant co-variates with increased data 14th June 2019

fPlotGAMMs(m1_larv, "Larvaceans")
dev.copy2pdf(file = "Larvaceans_sqrtsqrt.pdf", paper = "A4r")

saveRDS(m1_larv, file = "ModelOutput/gamm_Larvaceans.rds")
rm(m1_larv)
# write2pdf(summary(gmm1_larv), "LarvSummary.pdf",quiet = FALSE, title = "Larvaceans") # This will only work if you have tex installed

########## Euphausiids ##########
## There are a lot of zeros so I want to see what sampling methods have the problem
dat_Euph <- dat_Euph %>% mutate(Mins = case_when(Tow == "H" ~ Min$Min[Min$MinGrp == "Euphausiids_H"], 
                                                 Tow == "O" ~ Min$Min[Min$MinGrp == "Euphausiids_O"], 
                                                 Tow == "S1" ~ Min$Min[Min$MinGrp == "Euphausiids_S1"],
                                                 Tow == "S2" ~ Min$Min[Min$MinGrp == "Euphausiids_S2"],
                                                 Tow == "S3" ~ Min$Min[Min$MinGrp == "Euphausiids_S3"],
                                                 Tow == "V" ~ Min$Min[Min$MinGrp == "Euphausiids_V"])
)

zer_Euph2 <- dat_Euph %>% 
  group_by(Tow) %>% 
  summarise(PropZero = sum(TotAbundance==0) / n(), N = n())

zer_Euph2
# NAtl CPR is 3.1% euphausiids in data, but Richardson et al. (2006) says 87,236/207,619 = 42%
# So remove
dat_Euph <- dat_Euph %>% filter(Project!="SAHFOS" & Tow != "H") #
dat_Euph <- dat_Euph %>% filter(Tow != "V") #

# Lots of zeros drag the right hand down
ggplot(dat=dat_Euph, aes(x=log10(Chl), y=sqrt(TotAbundance))) + 
  facet_wrap(~ Tow, scales = "free", nrow = 1) +
  geom_point() +
  geom_smooth(method="lm")

ggplot(dat=dat_Euph, aes(x=log10(Chl), y=log10(TotAbundance + Mins))) + 
  facet_wrap(~ Tow, scales = "free", nrow = 1) +
  geom_point() +
  geom_smooth(method="lm")

# Plot with zeros removed
ggplot(dat=dat_Euph, aes(x=log10(Chl+0.0001), y=log10(TotAbundance+0.0001))) + 
  facet_wrap(~ Tow, scales = "free", nrow = 1) +
  geom_point() +
  geom_smooth(method="lm")

A <- dat_Euph %>% group_by(Tow, Project) %>% summarise(N = n(), Present = sum(TotAbundance > 0), Absent = sum(TotAbundance == 0), PropnPresent = Present / N)
# Look at these Projects vs Chl-a

ggplot(data = dat_Euph %>% filter(Tow == "H"), aes(x = log10(Chl), y = sqrt(TotAbundance))) +
  facet_wrap(~ Project, scales = "free") +
  geom_point() +
  geom_smooth(method="lm")

ggplot(data = dat_Euph %>% filter(Tow == "O"), aes(x = log10(Chl), y = sqrt(TotAbundance))) +
  facet_wrap(~ Project, scales = "free") +
  geom_point() +
  geom_smooth(method="lm")

ggplot(data = dat_Euph %>% filter(Tow == "V"), aes(x = log10(Chl), y = sqrt(TotAbundance))) +
  facet_wrap(~ Project, scales = "free") +
  geom_point() +
  geom_smooth(method="lm")

ggplot(data = dat_Euph %>% filter(Type == "CPR"), aes(x = log10(Chl), y = sqrt(TotAbundance))) +
  facet_wrap(~ Project, scales = "free") +
  geom_point() +
  geom_smooth(method="lm")

dat_Euph %>% filter(Project == "GILL Zooplankton Collection")
ggplot(data = dat_Euph %>% filter(Project == "GILL Zooplankton Collection"), aes(x = log10(Chl), y = sqrt(TotAbundance))) +
  facet_wrap(~ Transect, scales = "free") +
  geom_point() +
  geom_smooth(method="lm")
# Check why all negative for GILL data?

# CPR data dodgy - says 97
# The 3 options below are just some things Ant and I are playing around with to see if we can 
# Bring Chl and Bathy back into the models
X <- 10 # X% that CPR should be of analysis; (1-X) for Net - puts emphasis on Nets

# Check
CPR <- table(dat_Euph$Type)[1] * min(dat_Euph$WtVec)
Net <- table(dat_Euph$Type)[2] * max(dat_Euph$WtVec)
CPR
Net
CPR + Net

dat_Euph <- fAddWeighting2(dat_Euph, X)
CPR <- table(dat_Euph$Type)[1] * min(dat_Euph$WtVec2)
Net <- table(dat_Euph$Type)[2] * max(dat_Euph$WtVec2)
CPR
Net
CPR + Net


dat_Euph <- fAddWeighting(dat, "Euphausiids")
dat_Euph <- dat_Euph %>% filter(!(TotAbundance == 0 & Type == "CPR")) # All sig

m1_euph <- gam(log10(TotAbundance + MinEuph) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 fHarmonic(HarmDOY, k = 2) * fHarmonic(HarmHour, k = 2) + 
                 s(Longhurst, bs = "re") + s(Transect, bs = "re"),
               data = dat_Euph, weights = WtVec)


m1_euph <- gam(sqrt(sqrt(TotAbundance)) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                 # s(Longhurst, bs = "re") + s(Transect, bs = "re"),
                 s(Longhurst, bs = "re") + s(Project, bs = "re"),
               data = dat_Euph, weights = WtVec)

summary(m1_euph) # r2 = 53.5 %, log; r2 = 48.9% (with s(Project)), r2 = 56.1% (with s(Transect) sqrt(sqrt())
anova(m1_euph) # Chl and Bathy n.s 

m1_euph <- update(m1_euph, . ~ . -log10(Chl)) #
summary(m1_euph) # r2 = 54.9%
anova(m1_euph) # Bathy still n.s.

m1_euph <- update(m1_euph, . ~ . -s(Bathy, k = 3)) 
summary(m1_euph) # r2 = 54.9%
anova(m1_euph)  # All significant

fPlotGAMMs(m1_euph, "Euphausiids")
dev.copy2pdf(file = "Euphausiids_sqrtsqrt_Transect.pdf", paper = "A4r")


saveRDS(m1_euph, file = "ModelOutput/gamm_Euphausiids.rds")
rm(m1_euph)
# write2pdf(summary(gmm1_euph), "EuphSummary.pdf",quiet = FALSE, title = "Euphausiids") # This will only work if you have tex installed


########## Salps ########## 
A <- dat_Salp %>% group_by(Tow, Project) %>% summarise(N = n(), Present = sum(TotAbundance > 0), Absent = sum(TotAbundance == 0), PropnPresent = Present / N)
# Look at these Projects vs Chl-a

ggplot(data = dat_Salp %>% filter(Tow == "H" | Type == "CPR"), aes(x = log10(Chl), y = sqrt(TotAbundance))) +
  facet_wrap(~ Project, scales = "free") +
  geom_point() +
  geom_smooth(method="lm")

ggplot(data = dat_Salp %>% filter(Tow == "O"), aes(x = log10(Chl), y = sqrt(TotAbundance))) +
  facet_wrap(~ Project, scales = "free") +
  geom_point() +
  geom_smooth(method="lm")

ggplot(data = dat_Salp %>% filter(Tow == "V"), aes(x = log10(Chl), y = sqrt(TotAbundance))) +
  facet_wrap(~ Project, scales = "free") +
  geom_point() +
  geom_smooth(method="lm")


ggplot(data = dat_Salp %>% filter(Type == "CPR"), aes(x = log10(Chl), y = log10(TotAbundance+MinSalp))) +
  facet_wrap(~ Project, scales = "free") +
  geom_point() +
  geom_smooth(method="lm")

dat_Salp <- fAddWeighting2(dat_Salp, 0.1)

m1_salps <- gam(log10(TotAbundance + MinSalp) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) +
                  Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                  s(Longhurst, bs = "re") + s(Transect, bs = "re"),
                data = dat_Salp, weights = WtVec)
summary(m1_salps) # r2 = 55.9 %

m1_salps <- gam(sqrt(sqrt(TotAbundance)) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) +
                  Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                  s(Longhurst, bs = "re") + s(Transect, bs = "re"),
                data = dat_Salp, weights = WtVec2)
summary(m1_salps) # r2 = 49 %
anova(m1_salps) # all sig
fPlotGAMMs(m1_salps, "Salps")
dev.copy2pdf(file = "Salps_sqrtsqrt_Transect_CPR10percent.pdf", paper = "A4r")

visreg(m1_euph, "DOY2", by = "SST", breaks = c(0, 10, 20, 30), xlab = "Time of day", ylab = "log10(Abundance)", 
       scale = "response", overlay = TRUE, rug = 0, strip.names = c("0","10", "20", "30"))
visreg(m1_euph, "DOY2", by = "SST", breaks = c(0, 10, 20, 30), xlab = "Time of day", ylab = "log10(Abundance)", 
       scale = "response", overlay = TRUE, rug = 0, strip.names = c("0","10", "20", "30"))
visreg(m1_euph, "HarmHour", by = "HarmDOY", breaks = c(81/365*2*pi, 181/365*2*pi, 264/365*2*pi, 355/365*2*pi), xlab = "Time of Day", ylab = "log10(Abundance)", 
       scale = "response", overlay = TRUE, rug = 0, strip.names = c("Spr Equ","Sum Sol", "Aut Equ", "Win Sol"))


saveRDS(m1_salps, file = "ModelOutput/gamm_Salps.rds")
rm(m1_salps)
# write2pdf(summary(gmm1_salps), "SalpsSummary.pdf",quiet = FALSE, title = "Salps") # This will only work if you have tex installed


########## Jellyfish ########## 
# With zeros
ggplot(dat=dat_Jelly, aes(x=log10(Chl), y=log10(TotAbundance + MinJelly))) + 
  facet_wrap(~ Tow, scales = "free", nrow = 1) +
  geom_point() +
  geom_smooth(method="lm")

ggplot(dat=dat_Jelly, aes(x=log10(Chl), y=sqrt(TotAbundance))) + 
  facet_wrap(~ Tow, scales = "free", nrow = 1) +
  geom_point() +
  geom_smooth(method="lm")

ggplot(dat=dat_Euph, aes(x=log10(Chl), y=log10(TotAbundance + Mins))) + 
  facet_wrap(~ Tow, scales = "free", nrow = 1) +
  geom_point() +
  geom_smooth(method="lm")

# Plot with zeros removed
ggplot(dat=dat_Euph %>% filter(TotAbundance > 0), aes(x=log10(Chl+0.0001), y=log10(TotAbundance+0.0001))) + 
  facet_wrap(~ Tow, scales = "free", nrow = 1) +
  geom_point() +
  geom_smooth(method="lm")
m1_jelly <- gam(log10(TotAbundance + MinJelly) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) + 
                  Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                  s(Longhurst, bs = "re") + s(Transect, bs = "re"),
                data = dat_Jelly, weights = WtVec)

summary(m1_jelly) # r2 = 67.1 %
anova(m1_jelly) # All Sig

fPlotGAMMs(m1_jelly, "Jellyfish")

saveRDS(m1_jelly, file = "ModelOutput/gamm_Jellyfish.rds")


dat_Jelly2 <- dat_Jelly %>% filter(Type == "Net")
m2_jelly <- gam(log10(TotAbundance + MinJelly) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) + 
                  Mesh + Tow + exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 2) +
                  s(Longhurst, bs = "re") + s(Transect, bs = "re"),
                data = dat_Jelly2, weights = WtVec)
summary(m2_jelly) # r2 = 46.2
anova(m2_jelly) # All Sig

PlotGAMMs(m2_jelly, "Jellyfish_NoCPR")

rm(m1_jelly)
# write2pdf(summary(gmm1_jelly), "JellySummary.pdf",quiet = FALSE, title = "Jellyfish") # This will only work if you have tex installed

