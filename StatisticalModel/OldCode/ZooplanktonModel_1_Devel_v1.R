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


## PROBLEMS
# Check all the models, in parti


library(arsenal)
library(effects)
library(tidyverse)
library(lme4)
library(MuMIn)
library(splines)
library(mgcv)
library(visreg)


######################## Harmonic ##############################
source("~/Rcode/Harmonic.R")


dat <- readRDS("LatestDatabaseOuput_Final_Enviro.rds")

dat <- within(dat, {
  HarmHour <- (TimeLocal/24)*2*pi # Convert to radians
  HarmDOY <- (DOY2/365)*2*pi # Convert to radians
  Tow2 <- Tow
  Tow2[Tow =="S" | Tow == "H"] <- "H" 
  # Mesh[Mesh>570] <- 570
  Latitude2 <- abs(Latitude)
  TotAbundance[TotAbundance>10000] = 10000
})
dat <- droplevels(dat)

# For data already as frequencies
# ggplot(data = GenusFreq, aes(x = Var1, y = Freq)) +
#   geom_col() + theme_bw() +
#   scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
#                                    hjust = 1, face = "italic")) +
#   ylab("Frequency") + xlab("Genus")

# for data not as frequencies already
ggplot(data = dat, aes(x = log10(TotAbundance))) +
  geom_histogram() +
  theme_bw()







## These need to be moved to the databbse code:
dat$Gear_Mesh <- as.factor(dat$Gear_Mesh)

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

## GAM MODEL 1 WITH SPLINES FOR CONTINUOUS VARIABLES, 
##GEAR AND TOW TYPE INCLUDED AS FACTORS

# log10(Chl) + Harm(HarmHour,k=1) + exp(-Depth/1000) + Harm(HarmDOY,k=1)*poly(SST,2) + log10(Bathy) + Unit +
#   (1|User/Voyage) + offset(log(Vol)), data = dat)

m1 <- lm(log10(TotAbundance+min_val) ~ SST + log10(Chl) + log10(Bathy) +  Harm(HarmDOY,k=1)  + Gear_Mesh + exp(-Mid_Z/1000), 
         data = dat, subset = Group == "OmniCopepods")
summary(m1)
plot(allEffects(m1))

m2 <- lm(log10(TotAbundance+min_val) ~ poly(SST,2) + log10(Chl) + log10(Bathy) +  Harm(HarmDOY,k=1)  + Gear_Mesh + exp(-Mid_Z/1000), 
         data = dat, subset = Group == "OmniCopepods")
summary(m2)
plot(allEffects(m2))


m3 <- lm(log10(TotAbundance+min_val) ~ poly(SST,2) + log10(Chl) + log10(Bathy) + Gear + Mesh, 
         data = dat, subset = Group == "OmniCopepods")
summary(m3)
plot(allEffects(m3))

m4 <- lm(log10(TotAbundance+min_val) ~ ns(SST, df = 4) + log10(Chl) + log10(Bathy) + Harm(HarmDOY,k=1) + Mesh + Gear, 
         data = dat, subset = Group == "OmniCopepods")
summary(m4)
plot(allEffects(m4))



m5 <- lm(log10(TotAbundance+min_val) ~ ns(SST, df = 4) + log10(Chl) + log10(Bathy) + Harm(HarmDOY,k=1) + Gear + Mesh + Harm(HarmHour,k=1), 
         data = dat, subset = Group == "OmniCopepods")
summary(m5)
plot(allEffects(m5))


m6 <- lmer(log10(TotAbundance+min_val) ~ ns(SST, df = 4) + log10(Chl) + log10(Bathy) + Harm(HarmDOY,k=1) + Harm(HarmHour,k=1) + Mesh + (1|Gear), 
         data = dat, subset = Group == "OmniCopepods")
summary(m6)
plot(allEffects(m6))
r.squaredGLMM(m6)

ranef(m6,condVar=TRUE)

## Code for plotting random effects (From Max Carpenter)
REs <- ranef(m6, condVar = TRUE)
# Extract variance
qq <- attr(ranef(m6, condVar = TRUE)[[1]], "postVar")
# Extract intercepts
rand.interc <- REs$Gear
# Make a dataframe for plotting
df <- data.frame(Intercepts=REs$Gear[,1],
                 sd.interc=2*sqrt(qq[,,1:length(qq)]),
                 lev.names=rownames(rand.interc))
# Reorder levels
df$lev.names <-  factor(df$lev.names,rev(levels(df$lev.names)))
# 1200 * 900
ggplot(df,aes(lev.names,Intercepts)) + 
  geom_hline(yintercept=0, linetype = "dashed", color = "red") +
  geom_errorbar(aes(ymin=Intercepts-sd.interc, ymax=Intercepts+sd.interc),
                width = 0,color="black") +
  geom_point(aes(color = rev(lev.names)), size = 4) +
  guides(size = "none", shape = "none", color = "none") + theme_bw() +
  theme(axis.text.x=element_text(size=rel(1.2)), axis.title.x=element_text(size=rel(1.2)),
        axis.text.y=element_text(size=rel(1.4)), panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank()) +
  coord_flip() + ylab("Random intercept") + xlab("") +
  labs(title = "Random effect for Survey")





m7 <- lm(log10(TotAbundance+min_val) ~ ns(SST, df = 4) + log10(Chl) + log10(Bathy) 
         + Harm(HarmDOY,k=1) + Harm(HarmHour,k=1) + Mesh + Gear + Tow2, 
           data = dat, subset = Group == "OmniCopepods")
summary(m7)
plot(allEffects(m7))

# Didn't run
m8 <- lm(log10(TotAbundance+min_val) ~ ns(SST, df = 4) + log10(Chl) + log10(Bathy) 
         + Harm(HarmDOY,k=1) + Harm(HarmHour,k=1) + Mesh + Gear + Tow2 + ShipCruise, 
         data = dat, subset = Group == "OmniCopepods")
summary(m8)
plot(allEffects(m8))

# Bathy becomes insignificant
m9 <- lm(log10(TotAbundance+min_val) ~ ns(SST, df = 4) + log10(Chl) + log10(Bathy) 
         + Harm(HarmDOY,k=1) + Harm(HarmHour,k=1) + Mesh + Gear + Tow2 + Longhurst, 
         data = dat, subset = Group == "OmniCopepods")
summary(m9)
plot(allEffects(m9))

# 
m10 <- lm(log10(TotAbundance+min_val) ~ ns(SST, df = 4) + log10(Chl) + log10(Bathy) 
         + Harm(HarmDOY,k=1) + Harm(HarmHour,k=1) + Mesh + Gear + Tow2 + Longhurst + ns(Latitude, df = 2), 
         data = dat, subset = Group == "OmniCopepods")
summary(m10)
par(mfrow=c(2,2))
plot(allEffects(m10))
plot((m10))

m11 <- glm(TotAbundance+min_val ~ ns(SST, df = 4) + log10(Chl) + log10(Bathy) 
          + Harm(HarmDOY,k=1) + Harm(HarmHour,k=1) + Mesh + Gear + Tow2 + Longhurst 
          + ns(Latitude, df = 2), data = dat, subset = Group == "OmniCopepods", family=Gamma(link=log))
summary(m11)
plot(allEffects(m11))
plot(m11)



g1 <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=3,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
          Harm(HarmHour,k=1) + Mesh + Gear + Tow2 + Longhurst + ns(Latitude, df = 2),
          data = dat, subset = Group == "OmniCopepods")
summary(g1)

vis.gam(g1, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"

graphics.off()
quartz(width = 11, height = 8)

#png(filename="Net.png", units="in",  width=8, height=6, res=600)
par(mfrow = c(3,4))
par(mar = c(4,4,2,2))
visreg(g1, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(g1, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(g1, "HarmHour", scale = "response", xlab = "Time of Day")
visreg(g1, "Mesh", scale = "response", xlab = "Mesh  (microns)")
visreg(g1, "Gear", rug = FALSE, scale = "response", xlab = "Gear")
visreg(g1, "Tow2", rug = FALSE, scale = "response", xlab = "Tow")
visreg(g1, "Longhurst", rug = FALSE, scale = "response", xlab = "Longhurst")
visreg(g1, "Latitude", scale = "response", xlab = "Latitude")


# Surface for Time and Latitude - Not much there.....
g1a <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=3,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
            te(Latitude2,TimeLocal,k=3,bs=c("cr","cc")) + Mesh + Gear + Tow2 + Longhurst + ns(Latitude, df = 2),
          data = dat, subset = Group == "OmniCopepods")
summary(g1a)
plot(allEffects(g1a))
plot(g1a)
vis.gam(g1a, c("Latitude2", "TimeLocal"), type = "response", ticktype = "detailed", xlab = "\nLat", 
        ylab = "\nHour", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"

# Try with a spline, not a tensor
g2 <- gam(log10(TotAbundance+min_val) ~ s(SST,DOY2,k=20) + log10(Chl) + log10(Bathy) +
            Harm(HarmHour,k=1) + Mesh + Gear + Tow2 + Longhurst + ns(Latitude, df = 2),
          data = dat, subset = Group == "OmniCopepods")
summary(g2)
plot(allEffects(g2))
plot(g2)
vis.gam(g2, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"




## NEXT STEPS

# g1 seems to be the best model
# RANDOM EFFECTS - Use best model here and add Longhurst and Gear. Maybe ShipCruise
  # Gear as a random effect would be better as we don't need to choose a level in the maps
# HOW does the best model here (G1), work with the other functional groups

# gmm1 <- gamm(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=3,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
#             Harm(HarmHour,k=1) + Mesh + Gear + Tow2 + Longhurst + ns(Latitude, df = 2),
#            random=list(Gear=~1),
#           data = dat, subset = Group == "OmniCopepods")

########## Omnivorous Copepods ##########
gmm1_omni <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=3,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
            Harm(HarmHour,k=1) + Mesh + Tow2 + exp(-Mid_Z/1000) + Type +
              s(Gear, bs = "re") + s(Longhurst, bs = "re"),
          data = dat, subset = Group == "OmniCopepods")

# gmm1_omni <- lmer(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=3,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
#                    Harm(HarmHour,k=1) + Mesh + Tow2 + exp(-Mid_Z/1000) + Type +
#                    (1|Gear) + (1|Longhurst),
#                  data = dat, subset = Group == "OmniCopepods")






summary(gmm1_omni)

quartz(width = 11, height = 8)
vis.gam(gmm1_omni, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"
dev.print(pdf, 'Figures/OmniCope_Surface.pdf')

graphics.off()
quartz(width = 11, height = 8)
#png(filename="Net.png", units="in",  width=8, height=6, res=600)
par(mfrow = c(3,3))
par(mar = c(4,4,2,2))
visreg(gmm1_omni, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_omni, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_omni, "HarmHour", scale = "response", xlab = "Time of Day")
visreg(gmm1_omni, "Mesh", scale = "response", xlab = "Mesh  (microns)")
visreg(gmm1_omni, "Tow2", rug = FALSE, scale = "response", xlab = "Tow")
visreg(gmm1_omni, "Mid_Z", scale = "response", xlab = "Depth")
visreg(gmm1_omni, "Type", scale = "response", xlab = "Type")
########## 


dev.print(pdf, 'Figures/OmniCope.pdf')
saveRDS(gmm1_omni, file = "ModelOutput/gamm_Omni.rds")

write2pdf(summary(gmm1_omni), "OmniSummary.pdf",quiet = FALSE, title = "Omnivorous Copeods") # This will only work if you have tex installed

########## CarnCopepods ##########
# k = 4 needed for SST in this case
gmm1_carn <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=4,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                   Harm(HarmHour,k=1) + Mesh + Tow2 + exp(-Mid_Z/1000) + Type +
                   s(Gear, bs = "re") +  s(Longhurst, bs = "re"),
                 data = dat, subset = Group == "CarnCopepods")
summary(gmm1_carn)

graphics.off()
quartz(width = 11, height = 8)
vis.gam(gmm1_carn, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"
dev.print(pdf, 'Figures/CarnCope_Surface.pdf')

graphics.off()
quartz(width = 11, height = 8)
par(mfrow = c(3,3))
par(mar = c(4,4,2,2))
visreg(gmm1_carn, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_carn, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_carn, "HarmHour", scale = "response", xlab = "Time of Day")
visreg(gmm1_carn, "Mesh", scale = "response", xlab = "Mesh  (microns)")
visreg(gmm1_carn, "Tow2", rug = FALSE, scale = "response", xlab = "Tow")
# visreg(gmm1_carn, "Latitude", scale = "response", xlab = "Latitude")
visreg(gmm1_carn, "Mid_Z", scale = "response", xlab = "Depth")
visreg(gmm1_carn, "Type", scale = "response", xlab = "Type")

dev.print(pdf, 'Figures/CarnCope.pdf')

saveRDS(gmm1_carn, file = "ModelOutput/gamm_Carn.rds")

write2pdf(summary(gmm1_carn), "CarnSummary.pdf",quiet = FALSE, title = "Carnivorous Copeods") # This will only work if you have tex installed


########## Chaetognaths ##########
# k = 4 needed for SST in this case
gmm1_chaet <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=4,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                   Harm(HarmHour,k=1) + Mesh + Tow2 + exp(-Mid_Z/1000) + Type +
                   s(Gear, bs = "re") +  s(Longhurst, bs = "re"),
                 data = dat, subset = Group == "Chaetognaths")

summary(gmm1_chaet)
graphics.off()
quartz(width = 10, height = 7)
vis.gam(gmm1_chaet, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"

dev.print(pdf, 'Figures/Chaet_Surface.pdf')

quartz(width = 10, height = 7)
#png(filename="Net.png", units="in",  width=8, height=6, res=600)
par(mfrow = c(3,3))
par(mar = c(4,4,2,2))
visreg(gmm1_chaet, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_chaet, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_chaet, "HarmHour", scale = "response", xlab = "Time of Day")
visreg(gmm1_chaet, "Mesh", scale = "response", xlab = "Mesh  (microns)")
visreg(gmm1_chaet, "Tow2", rug = FALSE, scale = "response", xlab = "Tow")
visreg(gmm1_chaet, "Mid_Z", scale = "response", xlab = "Depth")
visreg(gmm1_chaet, "Type", scale = "response", xlab = "Type")

dev.print(pdf, 'Figures/Chaet.pdf')

saveRDS(gmm1_chaet, file = "ModelOutput/gamm_Chaet.rds")

write2pdf(summary(gmm1_chaet), "ChaetSummary.pdf",quiet = FALSE, title = "Chaetognaths") # This will only work if you have tex installed

########## Larvaceans ##########
# k = 4 needed for SST in this case
gmm1_larv <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=3,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                    Harm(HarmHour,k=1) + Mesh + Tow2 + exp(-Mid_Z/1000) + Type +
                    s(Gear, bs = "re") +  s(Longhurst, bs = "re"),
                  data = dat, subset = Group == "Larvaceans")

summary(gmm1_larv)
graphics.off()
quartz(width = 10, height = 7)
vis.gam(gmm1_larv, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"

dev.print(pdf, 'Figures/Larv_Surface.pdf')


quartz(width = 10, height = 7)
par(mfrow = c(3,3))
par(mar = c(4,4,2,2))
visreg(gmm1_larv, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_larv, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_larv, "HarmHour", scale = "response", xlab = "Time of Day")
visreg(gmm1_larv, "Mesh", scale = "response", xlab = "Mesh  (microns)")
visreg(gmm1_larv, "Tow2", rug = FALSE, scale = "response", xlab = "Tow")
visreg(gmm1_larv, "Mid_Z", scale = "response", xlab = "Depth")
visreg(gmm1_larv, "Type", scale = "response", xlab = "Type")

dev.print(pdf, 'Figures/Larv.pdf')
saveRDS(gmm1_larv, file = "ModelOutput/gamm_Larv.rds")

write2pdf(summary(gmm1_larv), "LarvSummary.pdf",quiet = FALSE, title = "Larvaceans") # This will only work if you have tex installed

########## Euphausiids ##########
# k = 4 needed for SST in this case
gmm1_euph <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=4,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                   Harm(HarmHour,k=1) + Mesh + Tow2 + exp(-Mid_Z/1000) + Type +
                   s(Gear, bs = "re") +  s(Longhurst, bs = "re"),
                 data = dat, subset = Group == "Euphausiids")

summary(gmm1_euph)
graphics.off()
quartz(width = 10, height = 7)
vis.gam(gmm1_euph, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"

dev.print(pdf, 'Figures/Euphaus_Surface.pdf')

quartz(width = 10, height = 7)
#png(filename="Net.png", units="in",  width=8, height=6, res=600)
par(mfrow = c(3,3))
par(mar = c(4,4,2,2))
visreg(gmm1_euph, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_euph, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_euph, "HarmHour", scale = "response", xlab = "Time of Day")
visreg(gmm1_euph, "Mesh", scale = "response", xlab = "Mesh  (microns)")
visreg(gmm1_euph, "Tow2", rug = FALSE, scale = "response", xlab = "Tow")
visreg(gmm1_euph, "Mid_Z", scale = "response", xlab = "Depth")
visreg(gmm1_euph, "Type", scale = "response", xlab = "Type")

dev.print(pdf, 'Figures/Euphaus.pdf')
saveRDS(gmm1_euph, file = "ModelOutput/gamm_Euphaus.rds")

write2pdf(summary(gmm1_euph), "EuphSummary.pdf",quiet = FALSE, title = "Euphausiids") # This will only work if you have tex installed

########## Jellyfish ########## 
# k = 4 needed for SST in this case
gmm1_jelly <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=4,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                   Harm(HarmHour,k=1) + Mesh + Tow2 + exp(-Mid_Z/1000) + Type +
                   s(Gear, bs = "re") +  s(Longhurst, bs = "re"),
                 data = dat, subset = Group == "Jellyfish")

summary(gmm1_jelly)
graphics.off()
quartz(width = 10, height = 7)
vis.gam(gmm1_jelly, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"

dev.print(pdf, 'Figures/Jelly_Surface.pdf')


quartz(width = 10, height = 7)
#png(filename="Net.png", units="in",  width=8, height=6, res=600)
par(mfrow = c(3,3))
par(mar = c(4,4,2,2))
visreg(gmm1_jelly, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_jelly, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_jelly, "HarmHour", scale = "response", xlab = "Time of Day")
visreg(gmm1_jelly, "Mesh", scale = "response", xlab = "Mesh  (microns)")
visreg(gmm1_jelly, "Tow2", rug = FALSE, scale = "response", xlab = "Tow")
visreg(gmm1_jelly, "Mid_Z", scale = "response", xlab = "Depth")
visreg(gmm1_jelly, "Type", scale = "response", xlab = "Type")

dev.print(pdf, 'Figures/Jelly.pdf')
saveRDS(gmm1_jelly, file = "ModelOutput/gamm_Jelly.rds")

write2pdf(summary(gmm1_jelly), "JellySummary.pdf",quiet = FALSE, title = "Jellyfish") # This will only work if you have tex installed

########## Salps ########## 
# k = 4 needed for SST in this case
gmm1_salps <- gam(log10(TotAbundance+min_val) ~ te(SST,DOY2,k=4,bs=c("cr","cc")) + log10(Chl) + log10(Bathy) +
                    Harm(HarmHour,k=1) + Mesh + Tow2 + exp(-Mid_Z/1000) + Type +
                    s(Gear, bs = "re") +  s(Longhurst, bs = "re"),
                  data = dat, subset = Group == "Salps")

summary(gmm1_salps)
graphics.off()
quartz(width = 10, height = 7)
vis.gam(gmm1_salps, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"

dev.print(pdf, 'Figures/Salps_Surface.pdf')
saveRDS(gmm1_salps, file = "ModelOutput/gamm_Salps.rds")

quartz(width = 10, height = 7)
#png(filename="Net.png", units="in",  width=8, height=6, res=600)
par(mfrow = c(3,3))
par(mar = c(4,4,2,2))
visreg(gmm1_salps, "Chl", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(gmm1_salps, "Bathy", scale = "response", xlab = "Bathy (m)")
visreg(gmm1_salps, "HarmHour", scale = "response", xlab = "Time of Day")
visreg(gmm1_salps, "Mesh", scale = "response", xlab = "Mesh  (microns)")
visreg(gmm1_salps, "Tow2", rug = FALSE, scale = "response", xlab = "Tow")
visreg(gmm1_salps, "Mid_Z", scale = "response", xlab = "Depth")
visreg(gmm1_salps, "Type", scale = "response", xlab = "Type")

dev.print(pdf, 'Figures/Salps.pdf')
saveRDS(gmm1_salps, file = "ModelOutput/gamm_Salps.rds")

write2pdf(summary(gmm1_salps), "SalpsSummary.pdf",quiet = FALSE, title = "Salps") # This will only work if you have tex installed

## Use the predict function to get the extreme values for chl for each group

# Dataframe for test prediction
test <- data.frame(SST = 15, DOY2 = 180, Chl = c(0.1, 5), Bathy = 1000, HarmHour = 1, Mesh = 100,
                   Tow2 = "H", Latitude = -30, Mid_Z = 0, Gear = 191, Longhurst = "ANTA", Type = "CPR")

omni <- predict(gmm1_omni,test)
carn <- predict(gmm1_carn,test)
chaet <- predict(gmm1_chaet,test)
euph <- predict(gmm1_euph,test)
salps <- predict(gmm1_salps,test)
jelly <- predict(gmm1_jelly,test)
larv <- predict(gmm1_larv,test)
# Try mapping the data and do the proportion maps of Ryans





graphics.off()
quartz(width = 11, height = 8)
#png(filename="Net.png", units="in",  width=8, height=6, res=600)
par(mfrow = c(3,4))
par(mar = c(4,4,2,2))
visreg(m1, "Types", rug = FALSE, scale = "response", xlab = "Type")
visreg(m1, "Net_ID_B", rug = FALSE, scale = "response", xlab = "Net ID")
visreg(m1, "Tow_ID", rug = FALSE, scale = "response", xlab = "Tow ID")
visreg(m1, "Mesh", scale = "response", xlab = "Mesh  (microns)")
visreg(m1, "Chla", scale = "response", xlab = "Chl-a (mg/m3)")
visreg(m1, "Depth", scale = "response", xlab = "Depth (m)")
visreg(m1, "Time2", scale = "response", xlab = "Time of Day")
visreg(m1, "SampleDepth", scale = "response", xlab = "Sample depth (m)")
visreg(m1, "SST", scale = "response", xlab = "SST")
visreg(m1, "DayOfYear", scale = "response", xlab = "Day of Year")

par(mar = c(1.5,0,1.5,0))
vis.gam(g1, c("SST", "DOY2"), type = "response", ticktype = "detailed", xlab = "\nSST (ºC)", 
        ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "cm", theta = 45, phi = 10, r = 100) # also "grey"









gm1 <- gam(log10(TOT_ABUND+min_val) ~ s(log10CHLO, fx = T, k = 5) +
             te(SST, day_of_year, bs = c("cr", "cc")) + s(BATHY, fx = T, k = 5) 
           + GEAR_MESH + s(TIMEloc), 
           data = krill_time)
summary(gm1)

library(effects)
plot(allEffects(m1))

## Save the GAM
gam_name = "krill_gam11.rds"
gam_save_file = paste(work_direct_main, "Zooplankton GAMs/Interpolate GAMs/", gam_name, sep = "")
saveRDS(gm1, file = gam_name) 
saveRDS(gm1, file = gam_save_file) # Save to Interpolate GAMs file

# PLOT RESIDUALS
par(mfrow = c(2,2))
gam.check(gm1)
abline(0,1, lwd = 2, col = "red")
par(mfrow = c(1,1))

## PLOT GAM
png(filename="krill_gam.png", units="in",  width=7, height=6, res=600)
par(mfrow = c(2,2), mar = c(3,3,0.5,0.5))
tt = visreg2d(gm1,  "SST",  "day_of_year", plot.type = "persp", theta = 50,
              col = colorRampPalette(brewer.pal(9,"YlOrRd"))(100), ylab = "\n Day of Year", xlab = "\n SST",
              zlab = "\n Predictor")
ano_pos = 0.17*max(tt$z)
mtext(text = "a)", las = 1, side = 2, line = 1.2, cex = 1.4, at = c(ano_pos), font = 1)


tt = visreg(gm1, xvar = "log10CHLO", ylab = "", yaxt = "n", xlab = "", xaxt = "n", partial = FALSE, rug = FALSE)
axis(1, mgp = c(0, 0.1, 0), cex.axis = 0.8, tck = -.01)
axis(2, mgp = c(0, 0.2, 0), cex.axis = 0.8, tck = -.01)
mtext(text = c(expression(paste("log"[10], "(Chlorophyll mg m"^-3, ")"))), side = 1, line = 1.2, cex = 0.9)
mtext(text = "Predictor", side = 2, line = 1, cex = 0.8, las = 3)
ano_pos = max(tt$fit[,"visregUpr"])
mtext(text = "b)", las = 2, side = 2, line = 1.2, cex = 1.4, at = c(ano_pos), font = 1)

tt= visreg(gm1, xvar = "BATHY", ylab = "", yaxt = "n", xlab = "",partial =FALSE, rug = FALSE)
#axis(1, mgp = c(0, 0.1, 0), cex.axis = 0.8, tck = -.01)
axis(2, mgp = c(0, 0.2, 0), cex.axis = 0.8, tck = -.01)
mtext(text = "Bathymetry (m)", side = 1, line = 1.2, cex = 0.9)
mtext(text = "Predictor", side = 2, line = 1.2, cex = 0.9, las = 3)
ano_pos = max(tt$fit[,"visregUpr"])
mtext(text = "c)", las = 2, side = 2, line = 1.2, cex = 1.4, at = c(ano_pos), font = 1)


tt=visreg(gm1, xvar = "GEAR_MESH", ylab = "", yaxt = "n", xlab = "", partial = FALSE, rug = FALSE)
num_gear = dim(tt$fit)[1]
#axis(1, labels = rep("",num_gear), at = seq(0.5*1/num_gear,1-0.5*1/num_gear,1/num_gear),mgp = c(0, 0.1, 0), cex.axis = 0.8, tck = -.01)
axis(2, mgp = c(0, 0.2, 0), cex.axis = 0.8, tck = -.01)
mtext(text = "Gear and Mesh Factor (GM)", side = 1, line = 1.2, cex = 0.9)
mtext(text = "Predictor", side = 2, line = 1.2, cex = 0.9, las = 3)
ano_pos = max(tt$fit[,"visregUpr"])
mtext(text = "d)", las = 2, side = 2, line = 1.2, cex = 1.4, at = c(ano_pos), font = 1)

dev.off()
