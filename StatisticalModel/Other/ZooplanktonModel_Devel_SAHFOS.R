# Rerunning the ZooplanktonModelDevel.R models so we can see how the CPR data looks.

library(arsenal)
library(effects)
library(tidyverse)
library(lme4)
library(MuMIn)
library(splines)
library(mgcv)
library(visreg)

dat <- readRDS("LatestDatabaseOuput_Final.rds")
dat$Mesh <- as.numeric(levels(dat$Mesh))[dat$Mesh]

dat <- within(dat, {
  HarmHour <- (TimeLocal/24)*2*pi # Convert to radians
  HarmDOY <- (DOY2/365)*2*pi # Convert to radians
  Tow2 <- Tow
  Tow2[Tow =="S" | Tow == "H"] <- "H" 
  Mesh[Mesh>570] <- 570
  Latitude2 <- abs(Latitude)
})

dat <- filter(dat,Project=="SAHFOS-CPR Atlantic Ocean")
dat <- droplevels(dat)


# We need a cruise ID but we don't have one in the current dataset. Temporarily create this so we can assess the 
# feasability of the random effects
dat <- dat %>% mutate(DeployID = paste0(ShipCruise, Year),
                      DeployID = paste0(DeployID, sprintf('%02g', Month)), # Now add in depth as well
                      DeployID = as.factor(DeployID)) 
 
######################## Harmonic ##############################
# Harmonic is to fit Hour and DOY
Harm <- function (theta, k = 4) {
  X <- matrix(0, length(theta), 2 * k)
  nam <- as.vector(outer(c("c", "s"), 1:k, paste, sep = ""))
  dimnames(X) <- list(names(theta), nam)
  m <- 0
  for (j in 1:k) {
    X[, (m <- m + 1)] <- cos(j * theta)
    X[, (m <- m + 1)] <- sin(j * theta)
  }
  X
}

############################  GAMs ############################  

# Models
# TotAbundance (Response, for each Group) ~ Type (fixed; CPR vs Net), Mesh (fixed linear), Tow (fixed; V, H, O, S), Gear (random), 
# ShipCruise (random), Longhurst (random), DOY (continuous), SST (continuous), Chl (continuous), 
# NewLocalTime (continuous), Bathy (continuous), maybe s(Latitude) + s(Longitude)
# Conditions: Mid_Z <= 200, all Years and Latitudes
# Error structure: Either gamma (log link function) or gaussian (log X + min_val)
# Use subset = ... in the glm/gam call for each Functional Group

# For log transformation
min_val = min(dat$TotAbundance[dat$TotAbundance > 0]) / 2

dat$DeployID = as.factor(dat$DeployID)
########## Omnivorous Copepods ##########
gmm1_omni <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=3,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                     Harm(HarmHour,k=1) + s(Longhurst, bs = "re") + s(DeployID, bs = "re"),
                   data = dat, subset = Group == "OmniCopepods")

summary(gmm1_omni)


graphics.off()
quartz(width = 11, height = 8)
par(mfrow = c(2,2))
par(mar = c(4,4,2,2))
visreg(gmm1_omni, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_omni, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_omni, "HarmHour", scale = "response", xlab = "Time of Day")
vis.gam(gmm1_omni, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"
########## 

dev.print(pdf, 'Figures/CPR_OmniCope.pdf')
saveRDS(gmm1_omni, file = "ModelOutput/CPR_gamm_Omni.rds")

write2pdf(summary(gmm1_omni), "CPR_OmniSummary.pdf",quiet = FALSE, title = "Omnivorous Copeods") # This will only work if you have tex installed

########## CarnCopepods ##########
# k = 4 needed for SST in this case
gmm1_carn <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=4,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                   Harm(HarmHour,k=1) + s(Longhurst, bs = "re") + s(DeployID, bs = "re"),
                 data = dat, subset = Group == "CarnCopepods")

summary(gmm1_carn)
graphics.off()
quartz(width = 11, height = 8)
par(mfrow = c(2,2))
par(mar = c(4,4,2,2))
visreg(gmm1_carn, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_carn, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_carn, "HarmHour", scale = "response", xlab = "Time of Day")
vis.gam(gmm1_carn, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"
dev.print(pdf, 'Figures/CPR_CarnCope.pdf')

saveRDS(gmm1_carn, file = "ModelOutput/CPR_gamm_Carn.rds")

write2pdf(summary(gmm1_carn), "CPR_CarnSummary.pdf",quiet = FALSE, title = "Carnivorous Copeods") # This will only work if you have tex installed

########## Chaetognaths ##########
# k = 4 needed for SST in this case
gmm1_chaet <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=4,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                   Harm(HarmHour,k=1) + s(Longhurst, bs = "re") + s(DeployID, bs = "re"),
                 data = dat, subset = Group == "Chaetognaths")

summary(gmm1_chaet)
graphics.off()

quartz(width = 10, height = 7)
par(mfrow = c(2,2))
par(mar = c(4,4,2,2))
visreg(gmm1_chaet, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_chaet, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_chaet, "HarmHour", scale = "response", xlab = "Time of Day")
vis.gam(gmm1_chaet, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"
dev.print(pdf, 'Figures/CPR_Chaet.pdf')

saveRDS(gmm1_chaet, file = "ModelOutput/CPR_gamm_Chaet.rds")

write2pdf(summary(gmm1_chaet), "CPR_ChaetSummary.pdf",quiet = FALSE, title = "Chaetognaths") # This will only work if you have tex installed

########## Larvaceans ##########
# k = 4 needed for SST in this case
gmm1_larv <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=3,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                    Harm(HarmHour,k=1) + s(Longhurst, bs = "re") + s(DeployID, bs = "re"),
                  data = dat, subset = Group == "Larvaceans")

summary(gmm1_larv)

graphics.off()

quartz(width = 10, height = 7)
par(mfrow = c(2,2))
par(mar = c(4,4,2,2))
visreg(gmm1_larv, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_larv, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_larv, "HarmHour", scale = "response", xlab = "Time of Day")
vis.gam(gmm1_larv, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"
dev.print(pdf, 'Figures/CPR_Larv.pdf')

saveRDS(gmm1_larv, file = "ModelOutput/CPR_gamm_Larv.rds")

write2pdf(summary(gmm1_larv), "CPR_LarvSummary.pdf",quiet = FALSE, title = "Larvaceans") # This will only work if you have tex installed

########## Euphausiids ##########
# k = 4 needed for SST in this case
gmm1_euph <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=4,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                   Harm(HarmHour,k=1) + s(Longhurst, bs = "re") + s(DeployID, bs = "re"),
                 data = dat, subset = Group == "Euphausiids")

summary(gmm1_euph)
  graphics.off()

quartz(width = 10, height = 7)
par(mfrow = c(2,2))
par(mar = c(4,4,2,2))
visreg(gmm1_euph, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_euph, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_euph, "HarmHour", scale = "response", xlab = "Time of Day")
vis.gam(gmm1_euph, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"

dev.print(pdf, 'Figures/CPR_Euphaus.pdf')
saveRDS(gmm1_euph, file = "ModelOutput/CPR_gamm_Euphaus.rds")

write2pdf(summary(gmm1_euph), "CPR_EuphSummary.pdf",quiet = FALSE, title = "Euphausiids") # This will only work if you have tex installed

########## Jellyfish ########## 
# k = 4 needed for SST in this case
gmm1_jelly <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=4,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                   Harm(HarmHour,k=1) + s(Longhurst, bs = "re") + s(DeployID, bs = "re"),
                 data = dat, subset = Group == "Jellyfish")

summary(gmm1_jelly)
graphics.off()

quartz(width = 10, height = 7)
par(mfrow = c(2,2))
par(mar = c(4,4,2,2))
visreg(gmm1_jelly, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_jelly, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_jelly, "HarmHour", scale = "response", xlab = "Time of Day")
vis.gam(gmm1_jelly, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"
dev.print(pdf, 'Figures/CPR_Jelly.pdf')

saveRDS(gmm1_jelly, file = "ModelOutput/CPR_gamm_Jelly.rds")

write2pdf(summary(gmm1_jelly), "CPR_JellySummary.pdf",quiet = FALSE, title = "Jellyfish") # This will only work if you have tex installed

########## Salps ########## 
# k = 4 needed for SST in this case
gmm1_salps <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=4,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                    Harm(HarmHour,k=1) + s(Longhurst, bs = "re") + s(DeployID, bs = "re"),
                  data = dat, subset = Group == "Salps")

summary(gmm1_salps)
graphics.off()

saveRDS(gmm1_salps, file = "ModelOutput/CPR_gamm_Salps.rds")

quartz(width = 10, height = 7)
par(mfrow = c(2,2))
par(mar = c(4,4,2,2))
visreg(gmm1_salps, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_salps, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_salps, "HarmHour", scale = "response", xlab = "Time of Day")
vis.gam(gmm1_salps, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"

dev.print(pdf, 'Figures/CPR_Salps.pdf')
saveRDS(gmm1_salps, file = "ModelOutput/CPR_gamm_Salps.rds")

write2pdf(summary(gmm1_salps), "CPR_SalpsSummary.pdf",quiet = FALSE, title = "Salps") # This will only work if you have tex installed

