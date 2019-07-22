# ZooplanktonModelDevel.R
# Ant, Jase, Ryan
# Created 11th April 2019
# Last Edited: 11th April 2019

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

library(effects)
library(tidyverse)
library(lme4)
library(MuMIn)
library(splines)
library(mgcv)
library(visreg)

source("fHarmonic.R") # Harmonic function

### FUNCTION ###
AddWeighting <- function (df, Grp) {
  # Takes a df and calculates weighting for each Type (Net vs CPR)
  Summary <- df %>% filter(Group == Grp) %>% group_by(Type) %>% 
    summarise(Num = n()) 
  Ratio <- Summary$Num[1] / Summary$Num[2]
  df <- df %>% filter(Group == Grp) %>% 
    mutate(WtVec = ifelse(Type == "Net", Ratio, 1), 
           WtVec = WtVec / mean(WtVec)) %>% 
    droplevels()
}
### FUNCTION ###

### FUNCTION ###
PlotGAMMs <- function (m1, Name) {
  Terms <- as.character(m1$terms)[3] # Terms from the model so we can print blank if n.s.
  
  x11(width = 15, height = 6)
  par(mfrow = c(2,5), mar = c(4,4,2,2))
  if(grepl("Chl", Terms, fixed = TRUE)) {
    visreg(m1, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")} else {
      plot.new()
      }
  if(grepl("Bathy", Terms, fixed = TRUE)) {
    visreg(m1, "Bathy", scale = "response", xlab = "Bathy (m)")} else {
      plot.new()
      }
  if(grepl("Mesh", Terms, fixed = TRUE)) {
    visreg(m1, "Mesh", scale = "response", xlab = "Mesh  (microns)")} else {
      plot.new()
      }
  
  if(grepl("Tow2", Terms, fixed = TRUE)) {
    visreg(m1, "Tow2", rug = FALSE, scale = "response", xlab = "Tow")} else {
      plot.new()
      }
  
  if(grepl("SST", Terms, fixed = TRUE)) {
    visreg(m1, "SST", scale = "response", xlab = "SST (ºC)")} else {
      plot.new()
      }
  
  if(grepl("DOY2", Terms, fixed = TRUE)) {
      visreg(m1, "DOY2", scale = "response", xlab = "Day of year")} else {
        plot.new()
        }
  
  # vis.gam(m1, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
  #        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100)
  visreg2d(m1, yvar = "DOY2", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
           ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)")
  if(grepl("Mid_Z", Terms, fixed = TRUE)) {
    visreg(m1, "Mid_Z", scale = "response", xlab = "Depth")} else {
      plot.new()
      }
  
  if(grepl("HarmHour", Terms, fixed = TRUE)) {
    visreg(m1, "HarmHour", scale = "response", xlab = "Time of Day")} else{
      plot.new()
      }
  
  # visreg(m1, "HarmHour", by = "Mid_Z", breaks = c(0, 100, 300, 1000), xlab = "Time of day", ylab = "log10(Abundance)", 
  #      scale = "response", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=100", "Depth=300", "Depth=1000"))
  
  if(grepl("exp(-Mid_Z/1000):Harm(HarmHour, k = 2)", Terms, fixed = TRUE)) {
    visreg(m1, "HarmHour", by = "Mid_Z", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", 
       scale = "response", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))} else {
         plot.new()
       }
  
  dev.print(pdf, paste0("Figures/", Name, ".pdf"))
}

# PlotGAMMs <- function (m1) {
# Can't get ggplot to work...
#   x11(width = 15, height = 6)  
#   p1 <- visreg(m1, "Chl", partial = FALSE, scale = "response", xlab = "Chl-a (mg/m3)", gg = TRUE) + theme_bw()
#   p2 <- visreg(m1, "Bathy", scale = "response", xlab = "Bathy (m)", gg = TRUE) + theme_bw()
#   p3 <- visreg(m1, "HarmHour", scale = "response", xlab = "Time of Day", gg = TRUE) + theme_bw()
#   p4 <- visreg(m1, "Mesh", scale = "response", xlab = "Mesh  (microns)", gg = TRUE) + theme_bw()
#   p5 <- visreg(m1, "Tow2", rug = FALSE, scale = "response", xlab = "Tow", gg = TRUE) + theme_bw()
#   p6 <- visreg(m1, "Mid_Z", scale = "response", xlab = "Depth", gg = TRUE) + theme_bw()
#   p7 <- visreg(m1, "SST", scale = "response", xlab = "SST (ºC)", gg = TRUE) + theme_bw()
#   p8 <- visreg(m1, "DOY2", scale = "response", xlab = "Day of year", gg = TRUE) + theme_bw()
#   # vis.gam(m1, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
#   #        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100)
#   # p9 <- visreg2d(m1_omni, yvar = "DOY2", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
#   #         ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", gg = TRUE) + theme_bw()
#   # p10 <- visreg(m1, "HarmHour", by = "Mid_Z", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", scale = "response", overlay = TRUE, gg = TRUE) + theme_bw()
#   library(gridExtra)
#   p <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8)
#   ggsave("Figures/OmniCope3.png", p, dpi = 1200)
#   # dev.print(pdf, 'Figures/OmniCope2.pdf')
# }

### FUNCTION ###

dat <- readRDS("LatestDatabaseOuput_Final_Enviro.rds") # n=977,881

dat <- within(dat, {
  HarmHour <- (TimeLocal/24)*2*pi # Convert to radians
  HarmDOY <- (DOY2/365)*2*pi # Convert to radians
  Tow2 <- Tow
  # Tow2[Tow =="S" | Tow == "H"] <- "H" 
  # Mesh[Mesh>570] <- 570
  Latitude2 <- abs(Latitude)
  TotAbundance[TotAbundance > 10000] = 10000
  Mid_Z[Mid_Z > 1000] = 1000
})
dat <- droplevels(dat)

### Look at data distribution
ggplot(data = dat, aes(x = log10(TotAbundance + 1))) +
  geom_histogram() +
  theme_bw()
summary(dat$TotAbundance) # No MD

ggplot(data = dat %>% filter(Type == "CPR"), aes(x = log10(TotAbundance+1))) +
  geom_histogram() +
  theme_bw()

ggplot(data = dat %>% filter(Type == "Net"), aes(x = log10(TotAbundance+1))) +
  geom_histogram() +
  theme_bw()

############################  GAMs ############################  
# RANDOM EFFECTS - Use best model here and add Longhurst and Gear. Maybe ShipCruise
# Gear as a random effect would be better as we don't need to choose a level in the maps
# HOW does the best model here (G1), work with the other functional groups

# Define min_val for GAM for each group - not worth doing it more elegantly!
Min <- dat %>% filter(TotAbundance > 0) %>% group_by(Group) %>% summarise(Min = min(TotAbundance) / 2)
MinChaet <- as.numeric(Min[1, 2])
MinCarn <- as.numeric(Min[2, 2])
MinOmni <- as.numeric(Min[3, 2])
MinEuph <- as.numeric(Min[4, 2])
MinJelly <- as.numeric(Min[5, 2])
MinLarv <- as.numeric(Min[6, 2])
MinSalp <- as.numeric(Min[7,2])

########## Omnivorous Copepods ##########
# For the tensor surface the different basis functions ("cr" = cubic regression, 
# "cs" = cubic regression with shrinkage, "tp" = thin plate spline, "ts" = tprs with shrinkage)
# all give the same results. Use "tp", as a bit better when gaps in data
# For s(Bathy), "tp" and "ts" have tighter SEs than "cr" abnd "cs"
# Bathy seems best with s(Bathy, k = 4)
# ti(SST, DOY2, k = 4, bs = c("tp","cc")) - k=4 appears best (less wiggly), but check
# Including Latitude with SST produces Latitude plot that doesn't really make sense (and gives very large SEs)
# s(SST, DOY2) gives a very folded surface - no good
# Harm(k=2) seems generally to be best
# NOTE: Harm(DOY2 * pi / 180, k = 2) * Harm(HarmHour, k = 2) is not useful
# Set up weighting vectors
# Each group has different weighting vector because each has different number of Net and CPR samples

dat_Omni <- AddWeighting(dat, "OmniCopepods")
dat_Carn <- AddWeighting(dat, "CarnCopepods")
dat_Chaet <- AddWeighting(dat, "Chaetognaths")
dat_Euph <- AddWeighting(dat, "Euphausiids")
dat_Jelly <- AddWeighting(dat, "Jellyfish")
dat_Larv <- AddWeighting(dat, "Larvaceans")
dat_Salp <- AddWeighting(dat, "Salps")

dat_Omni <- dat_Omni %>% mutate(Mid_Z = replace(Mid_Z, Mid_Z > 1000, 1000)) %>% droplevels()

m1_omni <- gam(log10(TotAbundance + MinOmni) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) + 
                 Mesh + Tow2 + exp(-Mid_Z/1000) * Harm(HarmHour, k = 2) +
                 # s(ShipCruise, bs = "re"),
                 s(Gear, bs = "re") + s(Longhurst, bs = "re") + s(Project, bs = "re"), 
               weights = WtVec, data = dat_Omni)
summary(m1_omni) #r2 = 53.6%
anova(m1_omni) # All significant
gam.check(m1_omni) # Looks reasonable I tink
PlotGAMMs(m1_omni, "CopepodOmnivores")
saveRDS(m1_omni, file = "ModelOutput/gamm_Omni.rds")
rm(m1_omni)

########## CarnCopepods ##########
m1_carn <- gam(log10(TotAbundance + MinCarn) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow2 + exp(-Mid_Z/1000) * Harm(HarmHour, k = 2) +
                 s(Gear, bs = "re") + s(Longhurst, bs = "re") + s(Project, bs = "re"), 
               weights = WtVec, data = dat_Carn)
summary(m1_carn) # r2 = 63.9%
anova(m1_carn) # Mesh n.s. so remove
m1_carn <- update(m1_carn, . ~ . -Mesh)
anova(m1_carn) # All significant (Tow2 = 0.0535)
summary(m1_carn) # r2 = 63.9%

PlotGAMMs(m1_carn, "CopepodCarnivores")
saveRDS(m1_carn, file = "ModelOutput/gamm_Carn.rds")
rm(m1_carn)
# write2pdf(summary(gmm1_carn), "CarnSummary.pdf",quiet = FALSE, title = "Carnivorous Copeods") # This will only work if you have tex installed

########## Chaetognaths ##########
# One chaetognath mesh size with 3000 um - make it next maximum mesh size of 710 for Chaets
dat_Chaet <- dat_Chaet %>% mutate(Mesh = replace(Mesh, Mesh > 710, 710))

m1_chaet <- gam(log10(TotAbundance + MinChaet) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) +
                  Mesh + Tow2 + exp(-Mid_Z/1000) * Harm(HarmHour, k = 2) +
                  s(Gear, bs = "re") +  s(Longhurst, bs = "re") + s(Project, bs = "re"),
                data = dat_Chaet, weights = WtVec)
summary(m1_chaet) # r2 = 52.1%
anova(m1_chaet) # Mesh n.s. so remove
m1_chaet <- update(m1_chaet, . ~ . -Mesh)
anova(m1_chaet) # All significant
summary(m1_chaet) # r2 = 52.1%
PlotGAMMs(m1_chaet, "Chaetognaths")
saveRDS(m1_chaet, file = "ModelOutput/gamm_Chaet.rds")
rm(m1_chaet)
# write2pdf(summary(gmm1_chaet), "ChaetSummary.pdf",quiet = FALSE, title = "Chaetognaths") # This will only work if you have tex installed

########## Larvaceans ##########
m1_larv <- gam(log10(TotAbundance + MinLarv) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow2 + exp(-Mid_Z/1000) * Harm(HarmHour, k = 2) +
                 s(Gear, bs = "re") +  s(Longhurst, bs = "re") + s(Project, bs = "re"),
               data = dat_Larv, weights = WtVec)
summary(m1_larv) # r2 = 35.4%
anova(m1_larv) # Interaction between Mid_Z and HarmHour n.s.
m1_larv <- update(m1_larv, . ~ . -exp(-Mid_Z/1000):Harm(HarmHour, k = 2))
anova(m1_larv) # Mesh n.s.
m1_larv <- update(m1_larv, . ~ . -Mesh)
summary(m1_larv) # r2 = 35.4%
PlotGAMMs(m1_larv, "Larvaceans")
saveRDS(m1_larv, file = "ModelOutput/gamm_Larv.rds")
rm(m1_larv)
# write2pdf(summary(gmm1_larv), "LarvSummary.pdf",quiet = FALSE, title = "Larvaceans") # This will only work if you have tex installed

########## Euphausiids ##########
# One euphausiid mesh size with 3000 um - make it next maximum mesh size of 710?
# Look at removing it???
dat_Euph <- dat_Euph %>% mutate(Mesh = replace(Mesh, Mesh > 710, 710))

m1_euph <- gam(log10(TotAbundance + MinEuph) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                 ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                 log10(Chl) + s(Bathy, k = 3) +
                 Mesh + Tow2 + exp(-Mid_Z/1000) * Harm(HarmHour, k = 2) +
                 s(Gear, bs = "re") +  s(Longhurst, bs = "re") + s(Project, bs = "re"),
               data = dat_Euph, weights = WtVec)
summary(m1_euph) # r2 = 44.8%
anova(m1_euph) # s(gear) and s(project) both very n.s.
m1_euph <- update(m1_euph, . ~ . -s(Gear, bs = "re")) # Then Project sig.
# m1_euph <- update(m1_euph, . ~ . -s(Project, bs = "re")) # Then Gear sig.
# ***Model has negative effect of Mesh when removing Gear (and negative n.s. Chl-a), 
# Model has positive effect of Mesh when removing Project (and negative sig. Chl-a) And removing 
anova(m1_euph) # Chl-a n.s.
m1_euph <- update(m1_euph, . ~ . -log10(Chl)) # Then Project sig.
anova(m1_euph) # All significant
summary(m1_euph) # r2 = 43.8%
PlotGAMMs(m1_euph, "Euphausiids")
saveRDS(m1_euph, file = "ModelOutput/gamm_Euphaus.rds")
rm(m1_euph)
# write2pdf(summary(gmm1_euph), "EuphSummary.pdf",quiet = FALSE, title = "Euphausiids") # This will only work if you have tex installed

########## Jellyfish ########## 
m1_jelly <- gam(log10(TotAbundance + MinJelly) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) + 
                  Mesh + Tow2 + exp(-Mid_Z/1000) * Harm(HarmHour, k = 2) +
                  s(Gear, bs = "re") +  s(Longhurst, bs = "re") + s(Project, bs = "re"),
                data = dat_Jelly, weights = WtVec)
summary(m1_jelly) # r2 = 62.1%
anova(m1_jelly) # Project n.s. so remove
m1_jelly <- update(m1_jelly, . ~ . -s(Project, bs = "re")) # Then Project sig.
anova(m1_jelly) # All sig.
summary(m1_jelly) # r2 = 61.8%
PlotGAMMs(m1_jelly, "Jellyfish")
saveRDS(m1_jelly, file = "ModelOutput/gamm_Jelly.rds")
rm(m1_jelly)
# write2pdf(summary(gmm1_jelly), "JellySummary.pdf",quiet = FALSE, title = "Jellyfish") # This will only work if you have tex installed

########## Salps ########## 
m1_salps <- gam(log10(TotAbundance + MinSalp) ~ # te(SST, DOY2, k = 4, bs = c("cr","cc")) + 
                  ti(SST, bs = "cr") + ti(DOY2, bs = "cc") + ti(SST, DOY2, k = 4, bs = c("cr","cc")) +
                  log10(Chl) + s(Bathy, k = 3) +
                  Mesh + Tow2 + exp(-Mid_Z/1000) * Harm(HarmHour, k = 2) +
                  s(Gear, bs = "re") +  s(Longhurst, bs = "re") + s(Project, bs = "re"),
                data = dat_Salp, weights = WtVec)
summary(m1_salps) # r2 = 50.9%
anova(m1_salps) # Everything significant
PlotGAMMs(m1_salps, "Salps")
saveRDS(m1_salps, file = "ModelOutput/gamm_Salps.rds")
rm(m1_salps)
# write2pdf(summary(gmm1_salps), "SalpsSummary.pdf",quiet = FALSE, title = "Salps") # This will only work if you have tex installed
