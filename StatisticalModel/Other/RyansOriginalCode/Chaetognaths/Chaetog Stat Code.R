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
library(rgdal)

work_direct_main = "~/Library/Mobile Documents/com~apple~CloudDocs/PhD/"# Set your main directory where you
                           # put "Zooplankton GAMs"

work_direct = paste(work_direct_main, "Zooplankton GAMs/Chaetognaths", sep = "")

setwd(work_direct)

chaetognaths <- within(read.csv("chaetognaths_2019.csv"),{
  
  TOT_ABUND = VALUE.per.volu
  
  LOWER_Z <- abs(LOWER_Z) # Make sure tow depth is positive
  UPPER_Z <- abs(UPPER_Z) # Make sure tow depth is positive
  
  MID_Z <- (UPPER_Z + LOWER_Z)/2 # Average tow depth
  MID_Z[MID_Z > 1000] <- 1000 # Max mean tow is 1000m (~30 obs greater than)
  
  MESH <- abs(MESH) # Make sure mesh size is positive
  BATHY <- abs(BATHY) # Make sure bathymetry is postive
  BATHY[BATHY > 5000] <- 5000 # Max depth is 5000 (~20 obs)
  CHLO[which(CHLO > 5)] <- 5 # Max chlo is 5 mg m^-3 (~120 obs)
  
  # Standardise Southern Hemisphere months to equivalent Northern Hemisphere months
  SH_1 <- c(LATITUDE < 0 & MON <=6) # Jan - June in South Hempishphere
  SH_2 <- (LATITUDE < 0 & MON  > 6) # July - December in South Hemisphere
  MON[SH_1] <- MON[SH_1] + 6
  MON[SH_2] <- MON[SH_2] - 6
  # 
  MON <- as.factor(MON) # Month as a factor
  
  ## DAY OF THE YEAR
  day_of_year <- round((as.integer(MON)-1)*30.4 + DAY)
  
  ## CLEAN UP TIME OF DAY
  TIMEloc[TIMEloc == -99.000] = NA
  TIMEloc <- TIMEloc %% 24
})


#source(paste(work_direct_main, "Zooplankton GAMs/Zoo Data Summary Plots.R", sep = ""))
#plot_zoo(chaetognaths, "chaetognaths_2019", save_plots = TRUE)

## REMOVE CPR DATA (~ 1000 obs)
#chaetognaths <- chaetognaths[-grep("CPR", chaetognaths$DATASET.ID),]
#chaetognaths <- chaetognaths[c(chaetognaths$GEAR != 191),]

## REMOVE ALL ORIG.VALUE THAT ISN'T # M^-3
chaetognaths <- chaetognaths[which(chaetognaths$Orig.UNITS %in% c("      #/m3", "    #/haul")),]

## REMOVE MESH SIZES OUTSIDE 100-500UM (~ 2000 obs)
#chaetognaths <- chaetognaths[c(chaetognaths$MESH >= 100 & chaetognaths$MESH <= 500),]
##

## REMOVE SAMPLES BEFORE 1958 (~ 4000 obs)
chaetognaths <- chaetognaths[c(chaetognaths$YEAR >= 1958),]

## Remove unecessary columns
chaetognaths <- chaetognaths[,-c(14:21,23,26:35,37:42)] # unnecessary columns
chaetognaths <- chaetognaths[!is.na(chaetognaths$TOT_ABUND),] #  presence/absence points

## CLEAN UP SPECIES NAMES
species_names <- as.character(chaetognaths$SCIENTIFIC.NAME....modifiers...)
species_namess <- sapply(strsplit(species_names, split = " -[", fixed = TRUE), function(x) (x[1]))
chaetognaths$SCIENTIFIC.NAME....modifiers...<- species_namess
colnames(chaetognaths)[17] <- "GROUP"

################# IMOS DATA ##########################
work_direct_imos = paste(work_direct_main, "Zooplankton GAMs/imos_2.csv", sep = "")
imos <- read.csv(work_direct_imos)
curr_zoo <- "Chaetognath"
wanted <- c("SHPCRUISE", "LATITUDE", "LONGITDE", "YEAR", "MON", "DAY", "TIMEloc", curr_zoo, "BATHY", "CHLO", "SST")
imos <- imos[,c(wanted)]
colnames(imos)[8] <- "TOT_ABUND"
imos$MESH <- 100
imos$GEAR <- 191 # CODED AS 'CPR' GEAR FOR ALL IMOS OBS - COPEPOD DATABASE GEAR NUMBER 191
imos$GROUP <- as.character("Chaetognatha")
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

chaetognaths$SHPCRUISE <- as.character(chaetognaths$SHPCRUISE)

chaetognaths <- merge(chaetognaths, imos, all = TRUE)

chaetognaths$SHPCRUISE <- as.factor(chaetognaths$SHPCRUISE)

## IDENTIFY WHICH ROWS ARE GENERAL "chaetognatha" VARIANTS
chaetognaths$CHAETOGNATHA <- FALSE

chaetognaths[grep("Chaetognatha", chaetognaths$GROUP), "CHAETOGNATHA"] <- TRUE

chaetognaths[is.na(chaetognaths$CHLO),"CHLO"] <- -50
chaetognaths[is.na(chaetognaths$SST),"SST"] <- -50
chaetognaths[is.na(chaetognaths$TIMEloc), "TIMEloc"] <- -50

chaetognaths <- chaetognaths[,-c(20:25)]

agg_chaetog <- aggregate(TOT_ABUND ~ SHPCRUISE + YEAR + MON + DAY
                       + LATITUDE + LONGITDE + T + MID_Z + 
                       + GEAR + MESH + CHAETOGNATHA + SST + CHLO + BATHY + MID_Z
                       + TIMEloc + day_of_year, 
                       data = chaetognaths, FUN = sum, na.action = NULL)

## For aggregated data where species are recorded as just "Chaetognatha", we remove 
## data from the same sample that has more specific taxa information. We do this 
## assuming that these more speciated observations would be included in the basic
## "Chaetognatha" samples, and we don't want to record them twice
#agg_chaetog <- agg_chaetog[order(agg_chaetog$CHAETOGNATHA, decreasing = TRUE),] ## PUT ALL CHAETOGNATHA OBSERVATIONS FIRST
#non_chaetognatha <- duplicated(agg_chaetog[,-c(13,19)]) ## IDENTIFY ALL CHAETOGNATH OBSERVATIONS FROM THE SAME
## SAMPLE AS CHAETOGNATHA OBSERVATIONS
#agg_chaetog <- agg_chaetog[c(!non_chaetognatha),] ## REMOVE ALL FALSE (NON-CHAETOGNATHA REPLICATES)

chaetognaths <- agg_chaetog

## PUT NA'S BACK
chaetognaths[chaetognaths$SST == -50, "SST"] <- NA
chaetognaths[chaetognaths$CHLO == -50, "CHLO"] <- NA
chaetognaths[chaetognaths$TIMEloc == -50, "TIMEloc"] <- NA

rm(agg_chaetog, non_chaetognatha)

chaetognaths$GEAR_MESH <- as.factor(paste(as.character(chaetognaths$GEAR), as.character(chaetognaths$MESH), sep = "_"))
chaetognaths$log10CHLO <- log10(chaetognaths$CHLO)

## ADD IN LONGHURST PROVINCES
setwd(paste(work_direct_main, "Zooplankton GAMs", sep = ""))
longhurst <- readOGR(dsn = "Longhurst", layer = "Longhurst_world_v4_2010")
curr_zoo <- chaetognaths
curr_locs <- data.frame("lat" = curr_zoo$LATITUDE, "lon" = curr_zoo$LONGITDE, 
                        "abund" = curr_zoo$TOT_ABUND)
coordinates(curr_locs) <- ~lon+lat
crs(curr_locs) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
provinces <- longhurst[,"ProvCode"]
a.data <- over(curr_locs, provinces)
curr_zoo$LONGHURST <- a.data$ProvCode
setwd(work_direct)

chaetognaths <- curr_zoo

#### CLEAN UP TIMES, SO THAT THEY ARE ALL NUMERIC
curr_zoo <- chaetognaths

fixtimes <- grep(":", curr_zoo$TIMEloc) # Which times are recorded as "HH:MM:SS"
fixtimes2 <- curr_zoo[fixtimes, "TIMEloc"] # Pull out these times

fixedtimes <- sapply(strsplit(fixtimes2,":"), # Convert to numeric equivalent
                     function(x) {
                       x <- as.numeric(x)
                       x[1]+x[2]/60
                     }
)

curr_zoo[fixtimes, 'TIMEloc'] <- as.character(fixedtimes)
chaetognaths <- curr_zoo

write.csv(chaetognaths, "clean_chaetogs.csv", row.names = FALSE)

## PLOT THE DATA - THIS TAKES A FEW MINUTES
#source(paste(work_direct_main, "Zooplankton GAMs/Zoo Data Summary Plots.R", sep = ""))
#plot_zoo(chaetognaths, "chaetognaths_2019c", save_plots = TRUE)

#####################################################################
#####################################################################

############## BUILD THE GAM

#####################################################################
#####################################################################

chaetognaths <- read.csv("clean_chaetogs.csv")
chaetog_time <- chaetognaths#[!is.na(chaetognaths$TIMEloc),] # If you only want obs with a local time stamp
chaetog_time[chaetog_time$MESH > 1000, 'MESH'] <- 1000 # Constrain mesh sizes to <= 1000um
chaetog_time <- chaetog_time[chaetog_time$MID_Z <= 200,] ## CHAETOGNATHS TOP 200 METRES (lose ~ 5% obs)

############################################################
## LAT CUT OFF - WHERE 15% OF THE DATA IS AFTER THIS LATITUDE
lat_cut = data.frame("lat" = seq(0, 90, 1), "num" = NA)
curr_zoo = chaetog_time
curr_zoo$latits = abs(curr_zoo$LATITUDE)

for(i in 1:dim(lat_cut)[1]){
  curr_lat = lat_cut[i,1]
  lat_cut[i,"num"] =  sum((curr_zoo$latits >= curr_lat), na.rm = TRUE) 
}

lat_cut = lat_cut[min(which(lat_cut[1:90,"num"] < 0.15*dim(curr_zoo)[1])), "lat"]

curr_zoo = curr_zoo[c(abs(curr_zoo$LATITUDE) < lat_cut),]
###################################################

chaetog_time = curr_zoo ## If you want to keep all the latitudes, don't run this line

min_val = min(chaetog_time[chaetog_time$TOT_ABUND > 0, 'TOT_ABUND'], na.rm =TRUE)/2

### MODEL WITH SPLINE SMOOTHERS FOR EVERYTHING
gm1 <- gam(log10(TOT_ABUND+min_val) ~  s(log10CHLO, fx = T, k = 3) +
             te(SST, day_of_year, bs = c("cr", "cc")) + s(BATHY, fx = T, k = 3) 
           + GEAR_MESH, 
           data = chaetog_time)
summary(gm1)

## SAVE GAM
gam_name = "chaets_gam12.rds"
gam_save_file = paste(work_direct_main, "Zooplankton GAMs/Interpolate GAMs/", gam_name, sep = "")
saveRDS(gm1, file = gam_name) 
saveRDS(gm1, file = gam_save_file) # Save to Interpolate GAMs file

# PLOT RESIDUALS
par(mfrow = c(2,2))
gam.check(gm1)
abline(0,1, lwd = 2, col = "red")
par(mfrow = c(1,1))

## PLOT THE MODEL
png(filename="chaetog_gam.png", units="in",  width=7, height=6, res=600)
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
