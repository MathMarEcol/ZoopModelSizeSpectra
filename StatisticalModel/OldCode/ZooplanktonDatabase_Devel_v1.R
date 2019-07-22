# ZooplanktonModels.R
# Ryan, Ant and Jase
# Last updated: 8/4/2019

library(mgcv)
library(effects)
library(splines)
library(visreg)
library(stringr)
library(tidyverse)
lapply(c("sp", "rgdal", "tmap"), library, character.only = TRUE)

dat1 <- read.csv("krill_2019.csv")
dat1 <- dat1 %>% mutate(Group = "Euphausiids")
dat2 <- read.csv("jellys_2019.csv")
dat2 <- dat2 %>% mutate(Group = "Jellyfish")
dat3 <- read.csv("chaetognaths_2019.csv")
dat3 <- dat3 %>% mutate(Group = "Chaetognaths")
dat4 <- read.csv("larvs_2019.csv")
dat4 <- dat4 %>% mutate(Group = "Larvaceans")
dat5 <- read.csv("salps_2019.csv")
dat5 <- dat5 %>% mutate(Group = "Salps")

# dat6 <- read.csv("copepods_2019.csv")
# dat7 <- read.csv("imos_2019.csv")
# dat8 <- read.csv("jedi_2019.csv")

dat <- rbind(dat1, dat2, dat3, dat4, dat5)
rm(dat1, dat2, dat3, dat4, dat5)

Retain <- c("Group", "SHPCRUISE", "YEAR", "MON", "DAY", "TIMEgmt", 
            "TIMEloc", "LATITUDE", "LONGITDE", "UPPER_Z", "LOWER_Z", 
            "T", "GEAR", "MESH", "VALUE.per.volu", "UNITS", 
            "SCIENTIFIC.NAME....modifiers...", "RECORD.ID", "DATASET.ID", "SHIP", "PROJ", "INST", "SST", "CHLO", "BATHY")
dat <- dat[, Retain]
dat <- dat %>% rename(Group = Group, ShipCruise = SHPCRUISE, Year = YEAR, Month = MON, 
                      Day = DAY, TimeGMT = TIMEgmt, TimeLoc = TIMEloc, 
                      Latitude = LATITUDE, Longitude = LONGITDE, Upper_Z = UPPER_Z, 
                      Lower_Z = LOWER_Z, TOW = T, Gear = GEAR, Mesh = MESH, 
                      Abundance = VALUE.per.volu, Units = UNITS, ScientificName = SCIENTIFIC.NAME....modifiers..., 
                      Record.ID = RECORD.ID, Dataset.ID = DATASET.ID, Ship = SHIP, Proj = PROJ, 
                      Inst = INST, SST = SST, Chl = CHLO, Bathy = BATHY)
dat <- dat %>% mutate(Diff = Upper_Z - Lower_Z)

# Check Bathy - Note: 0s must be where it is shallow - recode them to -20 m
summary(dat$Bathy)
sum(dat$Bathy < 0)
sum(dat$Bathy == 0)
sum(dat$Bathy > 0)
dat$Bathy[dat$Bathy == 0] <- -20
dat <- dat %>% mutate(Bathy = abs(Bathy), # Make Bathy positive
                      Bathy = replace(Bathy, Bathy > 7000, 7000)) # Set max bathymetry

# Check Lower_Z and Upper_Z
hist(dat$Diff)
sum(dat$Diff > 0)
# Remove samples where Lower_Z deeper than Upper_Z
hist(dat$Lower_Z)
dat <- dat %>% filter(Diff <= 0)
summary(dat$Lower_Z)
sum(dat$Lower_Z == -999.9)
hist(dat$Upper_Z)
summary(dat$Upper_Z)
sum(dat$Upper_Z == -999.9)
dat <- dat %>% filter(Lower_Z != -999.9) %>%
  mutate(Mid_Z = (Upper_Z + Lower_Z)/2) # Average tow depth
hist(dat$Mid_Z)

# Sort out mesh
hist(dat$Mesh)
summary(dat$Mesh)
sum(dat$Mesh == -999.0)
dat <- dat %>% filter(Mesh != -999.0) # Remove rows with no mesh size given

# Look at Chl
hist(dat$Chl)
dat <- dat %>% mutate(Chl = replace(Chl, Chl > 10, 10)) # Set max Chl

# Check Abundance
hist(log10(dat$Abundance+1)) # Looks OK

dat <- dat %>% mutate(UniqueSampleID = paste0(ShipCruise, Record.ID), # Create a Unique ID for every net deployment
                      UniqueSampleID = str_sub(UniqueSampleID, end = -6), # Remove the last 5 characters which are the 'taxaID')
                      UniqueSampleID = paste0(UniqueSampleID, sprintf('%04.0f', Mid_Z)), 
                      UniqueSampleID = as.factor(UniqueSampleID)) # Now add in the depth as well
dat <- dat %>% arrange(Group, UniqueSampleID)
dat2 <- dat %>% mutate(Group = as.factor(Group),
                       UniqueSampleID = as.factor(UniqueSampleID)) %>%
                       group_by(Group, UniqueSampleID) %>%
                       summarise(TotAbundance = sum(Abundance))
d <- dat %>%
  group_by(Group, UniqueSampleID) %>%
  distinct(Group, UniqueSampleID, .keep_all = TRUE)

dat3 <- right_join(d, dat2, by = c("Group", "UniqueSampleID"))







krill_uni <- summarize(krill, 
                       Abundance = sum(VALUE.per.volu,na.rm = TRUE),
                       n_combined = n(),
                       SHPCRUISE = SHPCRUISE[1],
                       YEAR = YEAR[1],
                       MON = MON[1],
                       DAY = DAY[1],
                       LATITUDE = LATITUDE[1],
                       LONGITDE = LONGITDE[1],
                       GEAR = GEAR[1],
                       MESH = MESH[1],
                       SST = SST[1],
                       CHLO = CHLO[1],
                       BATHY = BATHY[1],
                       MID_Z = MID_Z[1],
                       TIMEloc = TIMEloc[1],
                       day_of_year = day_of_year[1]
   
  # Standardise Southern Hemisphere months to equivalent Northern Hemisphere months
  SH_1 <- c(LATITUDE < 0 & MON <=6) # Jan - June in South Hempishphere
  SH_2 <- (LATITUDE < 0 & MON  > 6) # July - December in South Hemisphere
  MON[SH_1] <- MON[SH_1] + 6
  MON[SH_2] <- MON[SH_2] - 6

  MON <- as.factor(MON) # Month as a factor
  
  MID_Z <- (UPPER_Z + LOWER_Z)/2 # Average tow depth
  MID_Z[MID_Z > 1000] <- 1000 # Max mean tow is 1000m (~30 obs greater than)
  
  ## DAY OF THE YEAR
  day_of_year <- round((as.integer(MON)-1)*30.4 + DAY)
  
  ## CLEAN UP TIME OF DAY
  TIMEloc[TIMEloc == -99.000] = NA
  TIMEloc <- TIMEloc %% 24
})

#source(paste(work_direct_main, "Zooplankton GAMs/Zoo Data Summary Plots.R", sep = ""))
#plot_zoo(kt, "krill_2019", save_plots =TRUE)

## REMOVE ALL ORIG.VALUE THAT ISN'T # M^-3
krill <- krill[which(krill$Orig.UNITS %in% c("      #/m3", "    #/haul")),]

#krill <- krill[-grep("CPR", krill$DATASET.ID),]
#krill <- krill[c(krill$GEAR != 191),]

## Remove unecessary columns
krill <- krill[,-c(14:21,23,26:35,37:42)] # unnecessary columns

# Remove presence/absence points
krill <- krill[!is.na(krill$TOT_ABUND),] #  presence/absence points

## REMOVE SAMPLES FROM BEFORE 1958
# krill <- krill[c(krill$YEAR >= 1958),]

## CLEAN UP SPECIES NAMES
species_names <- as.character(krill$SCIENTIFIC.NAME....modifiers...)
species_namess <- sapply(strsplit(species_names, split = " -[", fixed = TRUE), function(x) (x[1]))
krill$SCIENTIFIC.NAME....modifiers...<- species_namess
colnames(krill)[17] <- "GROUP"

# Drop unused factor levels for SHPCRUISE
table(krill$SHPCRUISE)
nlevels(krill$SHPCRUISE)
krill <- droplevels(krill)


################# IMOS DATA ##########################
work_direct_imos = paste(work_direct_main, "Zooplankton GAMs/imos_2019.csv", sep = "") # OR "Zooplankton GAMs/imos_2.csv if you want to use 2017 imos data
imos <- read.csv(work_direct_imos)
table(imos$SHPCRUISE)

curr_zoo <- "Krill"
wanted <- c("PROJECT", "SHPCRUISE", "LATITUDE", "LONGITDE", "YEAR", "MON", "DAY", "TIMEloc", curr_zoo, "BATHY", "CHLO", "SST")
imos <- imos[,c(wanted)]
colnames(imos)[colnames(imos)=="Krill"] <- "TOT_ABUND"

imos$GEAR <- 191 # CODED AS 'CPR' GEAR FOR ALL IMOS OBS - COPEPOD DATABASE GEAR NUMBER 191
imos$GEAR[imos$PROJECT == "nrs"] <- 1000 # Code for NRS dropnet
imos$MESH <- 270
imos$MESH[imos$PROJECT == "nrs"] <- 100

imos$GROUP <- as.character("Euphausia")
imos$MID_Z <- 40 
imos$TIME <- as.factor("H")

library(lubridate)
imos$DATE <- ymd(sprintf('%04d%02d%02d',(imos$YEAR), imos$MON, imos$DAY, sep = ""),tz = "Europe/London")
imos$DOY <- yday(imos$DATE)

# Add different NRS to SHPCRUISE


# Standardise Southern Hemisphere months to equivalent Northern Hemisphere months
imos$SH_1 <- c(imos$LATITUDE < 0 & imos$MON <=6) # Jan - June in South Hempishphere
imos$SH_2 <- c(imos$LATITUDE < 0 & imos$MON  > 6) # July - December in South Hemisphere
imos$MON[imos$SH_1] <- imos$MON[imos$SH_1] + 6
imos$MON[imos$SH_2] <- imos$MON[imos$SH_2] - 6

hist(imos$BATHY)
imos$BATHY <- abs(imos$BATHY)

krill$SHPCRUISE <- as.character(krill$SHPCRUISE)
colnames(krill)[colnames(krill)=="day_of_year"] <- "DOY"

krill <- merge(krill, imos, all = TRUE)

krill$SHPCRUISE <- as.factor(krill$SHPCRUISE)

## IDENTIFY WHICH ROWS ARE GENERAL "EUPHAUSIA" VARIANTS
krill$KRILL_BASIC <- FALSE

krill[c(krill$GROUP == "Euphausia" | krill$GROUP ==  "Euphausia spp." | 
         krill$GROUP == "Euphausiacea"|krill$GROUP == "Euphausiacea spp." 
       | krill$GROUP == "Euphausiidae"),"KRILL_BASIC"] <- TRUE

krill[is.na(krill$CHLO),"CHLO"] <- -9999 # -9999 common MD code
krill[is.na(krill$SST),"SST"] <- -9999
krill[is.na(krill$TIMEloc), "TIMEloc"] <- -9999

# agg_krill <- aggregate(TOT_ABUND ~ SHPCRUISE + YEAR + MON + DAY
#                        + LATITUDE + LONGITDE + TIME +
#                          + GEAR + MESH + KRILL_BASIC + SST + CHLO + BATHY + MID_Z
#                        + TIMEloc + DOY, 
#                        data = krill, FUN = sum, na.action = NULL)

## For aggregated data where species are recorded as just "Euphausi" variant, we remove 
## data from the same sample that has more specific taxa information. We do this 
## assuming that these more speciated observations would be included in the basic
## "Euphausi" samples, and we don't want to record them twice
#agg_krill <- agg_krill[order(agg_krill$KRILL_BASIC, decreasing = TRUE),] ## PUT ALL KRILL_BASIC OBSERVATIONS FIRST
#non_krill_basic <- duplicated(agg_krill[,-c(13,19)]) ## IDENTIFY ALL NON-KRILL OBSERVATIONS FROM THE SAME
## SAMPLE AS COPEPODA AGGREGATE OBSERVATIONS
#agg_krill <- agg_krill[c(!non_krill_basic),] ## REMOVE ALL FALSE (NON-COPEPODA REPLICATES)

krill <- agg_krill

## PUT NA'S BACK
krill[krill$SST == -9999, "SST"] <- NA
krill[krill$CHLO == -9999, "CHLO"] <- NA
krill[krill$TIMEloc == -9999, "TIMEloc"] <- NA

rm(agg_krill)

## CREATE GEAR_MESH VARIABLE
krill$GEAR_MESH <- as.factor(paste(as.character(krill$GEAR), as.character(krill$MESH), sep =  "_"))

## LOG10CHLO VARIABLE
krill$log10CHLO <- log10(krill$CHLO)

## ADD IN LONGHURST PROVINCES
setwd(paste(work_direct_main, "Zooplankton GAMs", sep = ""))
longhurst <- readOGR(dsn = "Longhurst", layer = "Longhurst_world_v4_2010")
curr_zoo <- krill
curr_locs <- data.frame("lat" = curr_zoo$LATITUDE, "lon" = curr_zoo$LONGITDE, 
                         "abund" = curr_zoo$TOT_ABUND)
coordinates(curr_locs) <- ~lon+lat
crs(curr_locs) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
provinces <- longhurst[,"ProvCode"]
a.data <- over(curr_locs, provinces)
curr_zoo$LONGHURST <- a.data$ProvCode
setwd(work_direct)

krill <- curr_zoo

#### CLEAN UP TIMES, SO THAT THEY ARE ALL NUMERIC
curr_zoo <- krill

fixtimes <- grep(":", curr_zoo$TIMEloc) # Which times are recorded as "HH:MM:SS"
fixtimes2 <- curr_zoo[fixtimes, "TIMEloc"] # Pull out these times

fixedtimes <- sapply(strsplit(fixtimes2,":"), # Convert to numeric equivalent
                     function(x) {
                       x <- as.numeric(x)
                       x[1]+x[2]/60
                     }
)

curr_zoo[fixtimes, 'TIMEloc'] <- as.character(fixedtimes)
krill <- curr_zoo

write.csv(krill, "clean_krill.csv", row.names = FALSE)



###############################################################
#####################################################################

############## BUILD THE GAM

#####################################################################
#####################################################################

krill <- read.csv("clean_krill.csv")
krill_time <- krill#[!is.na(krill$TIMEloc),] ## only KRILL OBSERVATIONS WITH TIME STAMP
krill_time[c(krill_time$MESH > 1000), 'MESH'] <- 1000 ## limit mesh sizes to <= 1000um
krill_time <- krill_time[krill_time$MID_Z <= 200,] ## KRILL OBSERVATIONS FROM TOP 200M

############################################################
## LAT CUT OFF - WHERE 15% OF THE DATA IS AFTER THIS LATITUDE
lat_cut = data.frame("lat" = seq(0, 90, 1), "num" = NA)
curr_zoo = krill_time
curr_zoo$latits = abs(curr_zoo$LATITUDE)

for(i in 1:dim(lat_cut)[1]){
  curr_lat = lat_cut[i,1]
  lat_cut[i,"num"] =  sum((curr_zoo$latits >= curr_lat), na.rm = TRUE) 
}

lat_cut = lat_cut[min(which(lat_cut[1:90,"num"] < 0.15*dim(curr_zoo)[1])), "lat"]

curr_zoo = curr_zoo[c(abs(curr_zoo$LATITUDE) < lat_cut),]
################################################

krill_time = curr_zoo ## Don't run this line if you want all the lats

min_val = min(krill_time[krill_time$TOT_ABUND > 0, 'TOT_ABUND'], na.rm =TRUE)/2

## GAM MODEL 1 WITH SPLINES FOR CONTINUOUS VARIABLES, 
##GEAR AND TOW TYPE INCLUDED AS FACTORS
gm1 <- gam(log10(TOT_ABUND+min_val) ~ s(log10CHLO, fx = T, k = 5) +
             te(SST, day_of_year, bs = c("cr", "cc")) + s(BATHY, fx = T, k = 5) 
           + GEAR_MESH + s(TIMEloc), 
           data = krill_time)
summary(gm1)

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


# EXTRA
# %>%
#   filter(TotAbundance = sum(Abundance))
# 
# A <- dat %>% group_by(Group, UniqueSampleID) %>% distinct()
# 
# A <- distinct(dat, c(Group, UniqueSampleID), .keep_all = TRUE)
# 
# A <- dat %>% group_by(Group, UniqueSampleID) %>% filter(c == min(c))
# 
# 
# dat3 <- semi_join(dat2, dat, by = "UniqueSampleID")
# 
# A <- unique(inner_join(x = dat2, y = dat, by = "UniqueSampleID"), by = "UniqueSampleID")
# 
# dat3 <- base::merge(dat2, dat, by = "UniqueSampleID", all.y = TRUE)
# 
# 
# dat4 <- semi_join(dat, dat3, by = "UniqueSampleID" )
