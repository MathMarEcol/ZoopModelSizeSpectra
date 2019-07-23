rm(list=ls())
library(mgcv)
library(effects)
library(splines)
library(visreg)
library(raster)
library(maps)
library(mapdata)
library(maptools)
library(RColorBrewer)
library(stringr)
library(plyr)
library(RColorBrewer)
library(rgdal)

work_direct_main = "~/Library/Mobile Documents/com~apple~CloudDocs/PhD/"  # Set your main directory where you
# put "Zooplankton GAMs"

work_direct = paste(work_direct_main, "Zooplankton GAMs/Copepods", sep = "")

setwd(work_direct)

copepods <- within(read.csv("copepods_2019.csv"),{
  
  TOT_ABUND = VALUE.per.volu
  
  LOWER_Z <- abs(LOWER_Z) # Make sure tow depth is positive
  UPPER_Z <- abs(UPPER_Z) # Make sure tow depth is positive
  BATHY <- abs(BATHY) # Make sure bathymetry is postive
  BATHY[BATHY > 5000] <- 5000 # Max depth is 5000 (~20 obs)
  CHLO[which(CHLO > 5)] <- 5 # Max chlo is 5 mg m^-3 (~120 obs)
  
  # Standardise Southern Hemisphere months to equivalent Northern Hemisphere months
  SH_1 <- c(LATITUDE < 0 & MON <=6) # Jan - June in South Hempishphere
  SH_2 <- (LATITUDE < 0 & MON  > 6) # July - December in South Hemisphere
  MON[SH_1] <- MON[SH_1] + 6
  MON[SH_2] <- MON[SH_2] - 6
  # 
  
  MID_Z <- (UPPER_Z + LOWER_Z)/2 # Average tow depth
  MID_Z[MID_Z > 1000] <- 1000 # Max mean tow is 1000m (~30 obs greater than)
  
  MON <- as.factor(MON) # Month as a factor
  
  ## DAY OF THE YEAR
  day_of_year <- round((as.integer(MON)-1)*30.4 + DAY)
  
  ## CLEAN UP TIME OF DAY
  TIMEloc[TIMEloc == -99.000] = NA
  TIMEloc <- TIMEloc %% 24
  
})

### plot the data, this can take a few minutes
#source(paste(work_direct_main, "Zooplankton GAMs/Zoo Data Summary Plots.R", sep = ""))
#plot_zoo(copepods, "Copepods2019", save_plots = TRUE)

## CLEAN UP SPECIES NAMES
species_names <- as.character(copepods$SCIENTIFIC.NAME....modifiers...)
species_namess <- sapply(strsplit(species_names, split = " -[", fixed = TRUE), function(x) (x[1]))
copepods$SCIENTIFIC.NAME....modifiers...<- species_namess
colnames(copepods)[36] <- "GROUP"
copepods$GROUP <- as.factor(copepods$GROUP)

## REMOVE CPR DATA (~30,000 obs)
#copepods <- copepods[-grep("CPR", copepods$DATASET.ID),]
#copepods <- copepods[c(copepods$GEAR != 191),]

### REMOVE ALL ORIG.VALUE THAT ISN'T # M^-3
copepods <- copepods[which(copepods$Orig.UNITS %in% c("      #/m3", "    #/haul")),]

## REMOVE SAMPLES FROM BEFORE 1958 (~36,000 obs)
copepods <- copepods[c(copepods$YEAR >= 1958),]

## Remove unecessary columns 
copepods <- copepods[,-c(14:21,23,26:35,37:42)] # unnecessary columns

# Remove presence/absence points (~30,000 obs)
copepods <- copepods[!is.na(copepods$TOT_ABUND),]

## REMOVE MESH SIZES OUTSIDE 100-500UM (~25,000 obs)
#copepods <- copepods[c(copepods$MESH >= 100 & copepods$MESH <= 500),] #  mesh above 500um, below 0um (~1000 obs)
###

## IDENTIFY WHICH ROWS ARE GENERAL "COPEPODA" VARIANTS
copepods$COPEPODA <- FALSE
copepods[grep("Copepod", copepods$GROUP), "COPEPODA"] <- TRUE

################# IMOS DATA ##########################
## OMNI COPEPODS
work_direct_imos = paste(work_direct_main, "Zooplankton GAMs/imos_2019.csv", sep = "") # OR "Zooplankton GAMs/imos_2.csv if you want to use 2017 imos data
imos <- read.csv(work_direct_imos)
curr_zoo <- "Ominvore.Cope"
wanted <- c("SHPCRUISE", "LATITUDE", "LONGITDE", "YEAR", "MON", "DAY", "TIMEloc", curr_zoo, "BATHY", "CHLO", "SST")
imos <- imos[,c(wanted)]
colnames(imos)[8] <- "TOT_ABUND"
imos$MESH <- 100
imos$GEAR <- 191 # CODED AS 'CPR' GEAR FOR ALL IMOS OBS - COPEPOD DATABASE GEAR NUMBER 191
imos$GROUP <- as.character("Ominvore.Cope")
imos$MID_Z <- 100 
imos$T <- as.factor("H")
imos$day_of_year <- round((as.integer(imos$MON)-1)*30.4 + imos$DAY)
imos$SHPCRUISE <- as.character("IMOS")


# Standardise Southern Hemisphere months to equivalent Northern Hemisphere months
imos$SH_1 <- c(imos$LATITUDE < 0 & imos$MON <=6) # Jan - June in South Hempishphere
imos$SH_2 <- c(imos$LATITUDE < 0 & imos$MON  > 6) # July - December in South Hemisphere
imos$MON[imos$SH_1] <- imos$MON[imos$SH_1] + 6
imos$MON[imos$SH_2] <- imos$MON[imos$SH_2] - 6

imos$BATHY <- abs(imos$BATHY)
omni.imos <- imos

copepods$SHPCRUISE <- as.character(copepods$SHPCRUISE)

copepods <- merge(copepods, omni.imos, all = TRUE)

copepods$SHPCRUISE <- as.factor(copepods$SHPCRUISE)

## CARN COPEPODS
work_direct_imos = paste(work_direct_main, "Zooplankton GAMs/imos_2019.csv", sep = "") # OR "Zooplankton GAMs/imos_2.csv if you want to use 2017 imos data
imos <- read.csv(work_direct_imos)
curr_zoo <- "Carnivore.Cope"
wanted <- c("SHPCRUISE", "LATITUDE", "LONGITDE", "YEAR", "MON", "DAY", "TIMEloc", curr_zoo, "BATHY", "CHLO", "SST")
imos <- imos[,c(wanted)]
colnames(imos)[8] <- "TOT_ABUND"
imos$MESH <- 100
imos$GEAR <- 191 # CODED AS 'CPR' GEAR FOR ALL IMOS OBS - COPEPOD DATABASE GEAR NUMBER 191
imos$GROUP <- as.character("Carnivore.Cope")
imos$MID_Z <- 100 
imos$T <- as.factor("H")
imos$day_of_year <- round((as.integer(imos$MON)-1)*30.4 + imos$DAY)
imos$SHPCRUISE <- as.character("IMOS")


# Standardise Southern Hemisphere months to equivalent Northern Hemisphere months
imos$SH_1 <- c(imos$LATITUDE < 0 & imos$MON <=6) # Jan - June in South Hempishphere
imos$SH_2 <- c(imos$LATITUDE < 0 & imos$MON  > 6) # July - December in South Hemisphere
imos$MON[imos$SH_1] <- imos$MON[imos$SH_1] + 6
imos$MON[imos$SH_2] <- imos$MON[imos$SH_2] - 6

imos$BATHY <- abs(imos$BATHY)
carn.imos <- imos

copepods$SHPCRUISE <- as.character(copepods$SHPCRUISE)

copepods <- merge(copepods, carn.imos, all = TRUE)

copepods$SHPCRUISE <- as.factor(copepods$SHPCRUISE)

#### CLEAN UP TIMES, SO THAT THEY ARE ALL NUMERIC
curr_zoo <- copepods

fixtimes <- grep(":", curr_zoo$TIMEloc) # Which times are recorded as "HH:MM:SS"
fixtimes2 <- curr_zoo[fixtimes, "TIMEloc"] # Pull out these times

fixedtimes <- sapply(strsplit(fixtimes2,":"), # Convert to numeric equivalent
                     function(x) {
                       x <- as.numeric(x)
                       x[1]+x[2]/60
                     }
)

curr_zoo[fixtimes, 'TIMEloc'] <- as.character(fixedtimes)
copepods <- curr_zoo


######################################################################
### CREATE HERBIVORE/OMNIVORE/CARNIVORE COPEPOD DATA FRAME
feeding_type <- read.csv("Copepod_Feeding_Type.csv")

## REMOVE WHITESPACE FROM END OF NAMES
feeding_type[,c(1,2)] <- apply(feeding_type[,c(1,2)], 2, 
                               function(x){as.factor(str_trim(as.character(x), side = "right"))})

feeding_type[,1] <- as.factor(feeding_type[,1])

## REMOVE UNUSED COPEPOD SPECIES GROUPS
feeding_type[,1] <- factor(feeding_type[,1], levels(copepods$GROUP))
feeding_type[,2] <- as.factor(feeding_type[,2])
feeding_type <- feeding_type[complete.cases(feeding_type[,1]),]

copepods <- merge(copepods, feeding_type[,c("TAXON_NAME", "FEED")], 
                   by.x = "GROUP", by.y = "TAXON_NAME")


copepods$GROUP <-  mapvalues(copepods$GROUP, from = as.character(feeding_type[,1]),
                             to = as.character(feeding_type[,2]))

write.csv(copepods, "copepodss.csv", row.names = FALSE)

## AGGREGATE OMNIVORES AND HERBIVORES
oh_cops <- copepods[copepods$FEED %in% c("CH", "CO", "CBF", "CD", "CP", "SF"),]
  
### REMOVE NA'S FOR AGGREGATE FUNCTION
oh_cops[is.na(oh_cops$CHLO),"CHLO"] <- -50
oh_cops[is.na(oh_cops$SST),"SST"] <- -50
oh_cops[is.na(oh_cops$TIMEloc), "TIMEloc"] <- -50

### AGGREGATE OBSERVATIONS ACROSS SPECIES FROM SAME SAMPLE
agg_oh_cops <- aggregate(TOT_ABUND ~ SHPCRUISE + YEAR +MON + DAY+LATITUDE + LONGITDE+
                           T+MESH+GEAR+SST+CHLO+MID_Z+TIMEloc+BATHY +day_of_year, 
                      data = oh_cops, FUN = sum, na.action = NULL)

oh_cops <- agg_oh_cops

### PUT NA'S BACK IN
oh_cops[oh_cops$SST == -50, "SST"] <- NA
oh_cops[oh_cops$CHLO == -50, "CHLO"] <- NA
oh_cops[oh_cops$TIMEloc == -50, "TIMEloc"] <- NA

oh_cops$GEAR_MESH <- as.factor(paste(as.character(oh_cops$GEAR), as.character(oh_cops$MESH), sep = "_"))
oh_cops$log10CHLO <- log10(oh_cops$CHLO)

## ADD IN LONGHURST PROVINCES
setwd(paste(work_direct_main, "Zooplankton GAMs", sep = ""))
longhurst <- readOGR(dsn = "Longhurst", layer = "Longhurst_world_v4_2010")
curr_zoo <- oh_cops
curr_locs <- data.frame("lat" = curr_zoo$LATITUDE, "lon" = curr_zoo$LONGITDE, 
                        "abund" = curr_zoo$TOT_ABUND)
coordinates(curr_locs) <- ~lon+lat
crs(curr_locs) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
provinces <- longhurst[,"ProvCode"]
a.data <- over(curr_locs, provinces)
curr_zoo$LONGHURST <- a.data$ProvCode
setwd(work_direct)

oh_cops <- curr_zoo

write.csv(oh_cops, "clean_oh_cops.csv", row.names = FALSE)

## PLOT THE DATA - THIS TAKES A FEW MINUTES
#source(paste(work_direct_main, "Zooplankton GAMs/Zoo Data Summary Plots.R", sep = ""))
#plot_zoo(oh_cops, "omni cops", save_plots = FALSE)

## AGGREGATE CARNIVORES
carn_cops <- copepods[copepods$FEED %in% c("CC"),]

### REMOVE NA'S FOR AGGREGATE FUNCTION
carn_cops[is.na(carn_cops$CHLO),"CHLO"] <- -50
carn_cops[is.na(carn_cops$SST),"SST"] <- -50
carn_cops[is.na(carn_cops$TIMEloc), "TIMEloc"] <- -50

### AGGREGATE OBSERVATIONS ACROSS SPECIES FROM SAME SAMPLE
agg_carn_cops <- aggregate(TOT_ABUND ~ SHPCRUISE + YEAR + MON + DAY
                         + LATITUDE + LONGITDE  + T 
                         + GEAR + MESH + SST + CHLO + BATHY + MID_Z
                         + TIMEloc + day_of_year, 
                         data = carn_cops, FUN = sum, na.action = NULL)

carn_cops <- agg_carn_cops

### PUT NA'S BACK IN
carn_cops[carn_cops$SST == -50, "SST"] <- NA
carn_cops[carn_cops$CHLO == -50, "CHLO"] <- NA
carn_cops[carn_cops$TIMEloc == -50, "TIMEloc"] <- NA

carn_cops$GEAR_MESH <- as.factor(paste(as.character(carn_cops$GEAR), as.character(carn_cops$MESH), sep = "_"))
carn_cops$log10CHLO <- log10(carn_cops$CHLO)
carn_cops <- carn_cops[which(carn_cops$SST <= 30), ]

## ADD IN LONGHURST PROVINCES
setwd(paste(work_direct_main, "Zooplankton GAMs", sep = ""))
longhurst <- readOGR(dsn = "Longhurst", layer = "Longhurst_world_v4_2010")
curr_zoo <- carn_cops
curr_locs <- data.frame("lat" = curr_zoo$LATITUDE, "lon" = curr_zoo$LONGITDE, 
                        "abund" = curr_zoo$TOT_ABUND)
coordinates(curr_locs) <- ~lon+lat
crs(curr_locs) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
provinces <- longhurst[,"ProvCode"]
a.data <- over(curr_locs, provinces)
curr_zoo$LONGHURST <- a.data$ProvCode
setwd(work_direct)

carn_cops <- curr_zoo

write.csv(carn_cops, "clean_carn_cops.csv", row.names = FALSE)

######################################################################
######################################################################
### CREATE TOTAL_ABUNDANCE COPEPOD DATA FRAME ########################
######################################################################
######################################################################

copepods <- read.csv("copepodss.csv")

copepods$GEAR_MESH <- as.factor(paste(as.character(copepods$GEAR), as.character(copepods$MESH), sep = "_"))

copepods$TIMEloc <- as.character(copepods$TIMEloc)

### REMOVE NA'S FOR AGGREGATE FUNCTION
copepods[is.na(copepods$CHLO),"CHLO"] <- -50
copepods[is.na(copepods$SST),"SST"] <- -50
copepods[is.na(copepods$TIMEloc), "TIMEloc"] <- -50

copepods[is.na(copepods$COPEPODA), "COPEPODA"] <- FALSE

### AGGREGATE OBSERVATIONS ACROSS SPECIES FROM SAME SAMPLE
agg_cops <- aggregate(TOT_ABUND ~ SHPCRUISE + YEAR + MON + DAY
                      + LATITUDE + LONGITDE + T + MESH
                      + GEAR_MESH + COPEPODA + SST + CHLO + BATHY + MID_Z
                      + TIMEloc + day_of_year, 
                      data = copepods, FUN = sum, na.action = NULL)

## For aggregated data where species are recorded as just "Copepoda", we remove 
## data from the same sample that has more specific taxa information. We do this 
## assuming that these more speciated observations would be included in the basic
## "Copepoda" samples, and we don't want to record them twice

### REMOVE SPECIATED OBSERVATIONS FROM SAMPLE THAT HAS TOTAL COPEPOD SUM AND SPECIES BREAKDOWN
#agg_cops <- agg_cops[order(agg_cops$COPEPODA, decreasing = TRUE),] ## PUT ALL COPEPODA OBSERVATIONS FIRST

## IDENTIFY ALL NON-COPEPODA OBSERVATIONS FROM THE SAME SAMPLE AS COPEPODA AGGREGATE OBSERVATIONS
#non_copepoda <- duplicated(agg_cops[,-c(13,19)]) 

## REMOVE ALL FALSE (NON-COPEPODA REPLICATES)
#agg_cops <- agg_cops[c(!non_copepoda),]

copepods <- agg_cops

## PUT NA'S BACK
copepods[copepods$SST == -50, "SST"] <- NA
copepods[copepods$CHLO == -50, "CHLO"] <- NA
copepods[copepods$TIMEloc == -50, "TIMEloc"] <- NA


copepods$log10CHLO <- log10(copepods$CHLO)

## ADD IN LONGHURST PROVINCES
setwd(paste(work_direct_main, "Zooplankton GAMs", sep = ""))
longhurst <- readOGR(dsn = "Longhurst", layer = "Longhurst_world_v4_2010")
curr_zoo <- copepods
curr_locs <- data.frame("lat" = curr_zoo$LATITUDE, "lon" = curr_zoo$LONGITDE, 
                        "abund" = curr_zoo$TOT_ABUND)
coordinates(curr_locs) <- ~lon+lat
crs(curr_locs) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
provinces <- longhurst[,"ProvCode"]
a.data <- over(curr_locs, provinces)
curr_zoo$LONGHURST <- a.data$ProvCode
setwd(work_direct)

copepods <- curr_zoo

write.csv(copepods, "clean_copepod_tot_abundances.csv", row.names = FALSE)

## PLOT THE DATA - THIS TAKES A FEW MINUTES
#source(paste(work_direct_main, "Zooplankton GAMs/Zoo Data Summary Plots.R", sep = ""))
#plot_zoo(copepods, "total copepods", save_plots = FALSE)

#####################################################################
#####################################################################
#####################################################################

#####################################################################
#####################################################################

############## BUILD THE GAM (COPEPOD OMNIVORE/HERBIVORE MODEL)

#####################################################################
#####################################################################

omni.herb_cops <- read.csv("clean_oh_cops.csv")
oh_cop_time <- omni.herb_cops#[!is.na(omni.herb_cops$TIMEloc),] ## REMOVE NON LOCAL TIME-STAMP OBSERVATIONS 
oh_cop_time <- oh_cop_time[oh_cop_time$MID_Z <= 200,] ## TOP 200 METRES SAMPLES only 
oh_cop_time[oh_cop_time$MESH >1000, 'MESH'] <- 1000 # Constrain mesh sizes to <= 1000um

############################################################
## LAT CUT OFF - WHERE 15% OF THE DATA IS AFTER THIS LATITUDE
lat_cut = data.frame("lat" = seq(0, 90, 1), "num" = NA)
curr_zoo = oh_cop_time
curr_zoo$latits = abs(curr_zoo$LATITUDE)

for(i in 1:dim(lat_cut)[1]){
  curr_lat = lat_cut[i,1]
  lat_cut[i,"num"] =  sum((curr_zoo$latits >= curr_lat), na.rm = TRUE) 
}

lat_cut = lat_cut[min(which(lat_cut[1:90,"num"] < 0.15*dim(curr_zoo)[1])), "lat"]

curr_zoo = curr_zoo[c(abs(curr_zoo$LATITUDE) < lat_cut),]
###################################################

oh_cop_time = curr_zoo  ## If you want to keep all the latitudes, don't run this line

min_val = min(oh_cop_time[oh_cop_time$TOT_ABUND > 0, 'TOT_ABUND'], na.rm =TRUE)/2

### MODEL 1 WITH SPLINE SMOOTHERS FOR EVERYTHING
gm1 <- gam(log10(TOT_ABUND+min_val) ~  s(log10CHLO, fx = T, k = 5) +
             te(SST, day_of_year, bs = c("cr", "cc")) + s(BATHY, fx = T, k = 5) 
           + GEAR_MESH, 
           data = oh_cop_time) 
summary(gm1)

### PLOT RESIDUALS
par(mfrow = c(2,2))
gam.check(gm1)
abline(0,1, lwd = 2, col = "red")
par(mfrow = c(1,1))

### SAVE GAM
gam_name = "omnicops_gam12.rds"
gam_save_file = paste(work_direct_main, "Zooplankton GAMs/Interpolate GAMs/", gam_name, sep = "")
saveRDS(gm1, file = gam_name) 
saveRDS(gm1, file = gam_save_file) # Save to Interpolate GAMs file

## SAVE MODEL PLOTS
png(filename="ocop_gam.png", units="in",  width=7, height=6, res=600)
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




#####################################################################
#####################################################################
## COPEPOD CARNIVORE MODEL
#####################################################################
#####################################################################

carn_cops <- read.csv("clean_carn_cops.csv")
c_cop_time <- carn_cops#[!is.na(carn_cops$TIMEloc),] ## REMOVE NON LOCAL TIME-STAMP OBSERVATIONS
c_cop_time <- c_cop_time[c_cop_time$MID_Z <= 200,] ## TOP 200 METRES OBSERVATIONS only 
c_cop_time[c_cop_time$MESH >1000, 'MESH'] <- 1000 # Constrain mesh sizes to <= 1000um

############################################################
## LAT CUT OFF - WHERE 15% OF THE DATA IS AFTER THIS LATITUDE
lat_cut = data.frame("lat" = seq(0, 90, 1), "num" = NA)
curr_zoo = c_cop_time
curr_zoo$latits = abs(curr_zoo$LATITUDE)

for(i in 1:dim(lat_cut)[1]){
  curr_lat = lat_cut[i,1]
  lat_cut[i,"num"] =  sum((curr_zoo$latits >= curr_lat), na.rm = TRUE) 
}

lat_cut = lat_cut[min(which(lat_cut[1:90,"num"] < 0.15*dim(curr_zoo)[1])), "lat"]

curr_zoo = curr_zoo[c(abs(curr_zoo$LATITUDE) < lat_cut),]
###################################################

c_cop_time = curr_zoo ## If you want to keep all the latitudes, don't run this line

min_val = min(c_cop_time[c_cop_time$TOT_ABUND > 0, 'TOT_ABUND'], na.rm =TRUE)/2

### MODEL 1 WITH SPLINE SMOOTHERS FOR EVERYTHING
gm1 <- gam(log10(TOT_ABUND+min_val) ~  s(log10CHLO, fx = T, k =5) +
             te(SST, day_of_year, bs = c("cr", "cc")) + s(BATHY, fx = T, k = 5) 
           + GEAR_MESH, 
           data = c_cop_time)
summary(gm1)

## PLOT RESIDUALS
par(mfrow = c(2,2))
gam.check(gm1)
abline(0,1, lwd = 2, col = "red")
par(mfrow = c(1,1))

## SAVE GAM
gam_name = "copecarn_gam12.rds"
gam_save_file = paste(work_direct_main, "Zooplankton GAMs/Interpolate GAMs/", gam_name, sep = "")
saveRDS(gm1, file = gam_name) 
saveRDS(gm1, file = gam_save_file) # Save to Interpolate GAMs file

## SAVE MODEL PLOTS
png(filename="ccop_gam.png", units="in",  width=7, height=6, res=600)
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

#####################################################################
#####################################################################
#### MODEL FOR TOTAL ABUNDANCES OF COPEPODS (AGGREGATED DATASET)
#####################################################################
#####################################################################

copepods_total <- read.csv("clean_copepod_tot_abundances.csv")
copepods_total_t <- copepods_total#[!is.na(copepods_total$TIMEloc),]  ## REMOVE NON TIME STAMP DATA (~)
copepods_total_t <- copepods_total_t[copepods_total_t$MID_Z <= 200,] ## REMOVE OBS FROM DEEPER THAN 200 METRE SAMPLE DEPTH
copepods_total_t[copepods_total_t$MESH >1000, 'MESH'] <- 1000 ## MESH SIZE CONSTRAINED TO <= 1000 UM 

############################################################
## LAT CUT OFF - WHERE 15% OF THE DATA IS AFTER THIS LATITUDE
lat_cut = data.frame("lat" = seq(0, 90, 1), "num" = NA)
curr_zoo = copepods_total_t
curr_zoo$latits = abs(curr_zoo$LATITUDE)

for(i in 1:dim(lat_cut)[1]){
  curr_lat = lat_cut[i,1]
  lat_cut[i,"num"] =  sum((curr_zoo$latits >= curr_lat), na.rm = TRUE) 
}

lat_cut = lat_cut[min(which(lat_cut[1:90,"num"] < 0.1*dim(curr_zoo)[1])), "lat"]

curr_zoo = curr_zoo[c(abs(curr_zoo$LATITUDE) < lat_cut),]
###################################################

copepods_total_t = curr_zoo # Don't run this line if you want all the lats

min_val = min(copepods_total_t[copepods_total_t$TOT_ABUND > 0, 'TOT_ABUND'], na.rm =TRUE)/2

gm1 <- gam(log10(TOT_ABUND+min_val) ~  s(log10CHLO, fx = T, k = 5) +
             te(SST, day_of_year, bs = c("cr", "cc")) + s(BATHY, fx = T, k = 5) 
           + GEAR_MESH, data = copepods_total_t)
             
summary(gm1)

## PLOT RESIDUALS
par(mfrow = c(2,2))
gam.check(gm1)
abline(0,1, lwd = 2, col = "red")
par(mfrow = c(1,1))

## SAVE GAM
gam_name = "copetotal_gam12.rds"
gam_save_file = paste(work_direct_main, "Zooplankton GAMs/Interpolate GAMs/", gam_name, sep = "")
saveRDS(gm1, file = gam_name) 
saveRDS(gm1, file = gam_save_file) # Save to Interpolate GAMs file

### PLOT AND SAVE GAM
png(filename="totalcop_gam.png", units="in",  width=7, height=6, res=600)
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
