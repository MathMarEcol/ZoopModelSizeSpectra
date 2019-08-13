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
lapply(c("sp", "rgdal", "tmap"), library, character.only = TRUE)

work_direct_main ="~/Library/Mobile Documents/com~apple~CloudDocs/PhD/" # Set your main directory where you
# put "Zooplankton GAMs"

work_direct = paste(work_direct_main, "Zooplankton GAMs/Krill", sep = "")

setwd(work_direct)

krill <- within(read.csv("krill_2019.csv"),{
  
  TOT_ABUND = VALUE.per.volu
  
  LOWER_Z <- abs(LOWER_Z) # Make sure tow depth is positive
  UPPER_Z <- abs(UPPER_Z) # Make sure tow depth is positive
  MESH <- abs(MESH) # Make sure mesh size is positive
  BATHY <- abs(BATHY) # Make sure bathymetry is postive
  BATHY[BATHY > 5000] <- 5000 # Max depth is 5000 (~20 obs)
  CHLO[which(CHLO > 5)] <- 5 # Max chlo is 5 mg m^-3 (~120 obs)
  
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
krill <- krill[c(krill$YEAR >= 1958),]

## CLEAN UP SPECIES NAMES
species_names <- as.character(krill$SCIENTIFIC.NAME....modifiers...)
species_namess <- sapply(strsplit(species_names, split = " -[", fixed = TRUE), function(x) (x[1]))
krill$SCIENTIFIC.NAME....modifiers...<- species_namess
colnames(krill)[17] <- "GROUP"

################# IMOS DATA ##########################
work_direct_imos = paste(work_direct_main, "Zooplankton GAMs/imos_2019.csv", sep = "") # OR "Zooplankton GAMs/imos_2.csv if you want to use 2017 imos data
imos <- read.csv(work_direct_imos)
curr_zoo <- "Krill"
wanted <- c("SHPCRUISE", "LATITUDE", "LONGITDE", "YEAR", "MON", "DAY", "TIMEloc", curr_zoo, "BATHY", "CHLO", "SST")
imos <- imos[,c(wanted)]
colnames(imos)[8] <- "TOT_ABUND"
imos$MESH <- 100
imos$GEAR <- 191 # CODED AS 'CPR' GEAR FOR ALL IMOS OBS - COPEPOD DATABASE GEAR NUMBER 191
imos$GROUP <- as.character("Euphausia")
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

krill$SHPCRUISE <- as.character(krill$SHPCRUISE)

krill <- merge(krill, imos, all = TRUE)

krill$SHPCRUISE <- as.factor(krill$SHPCRUISE)

## IDENTIFY WHICH ROWS ARE GENERAL "EUPHAUSIA" VARIANTS
krill$KRILL_BASIC <- FALSE

krill[c(krill$GROUP == "Euphausia" | krill$GROUP ==  "Euphausia spp." | 
         krill$GROUP == "Euphausiacea"|krill$GROUP == "Euphausiacea spp." 
       | krill$GROUP == "Euphausiidae"),"KRILL_BASIC"] <- TRUE

krill[is.na(krill$CHLO),"CHLO"] <- -50
krill[is.na(krill$SST),"SST"] <- -50
krill[is.na(krill$TIMEloc), "TIMEloc"] <- -50

agg_krill <- aggregate(TOT_ABUND ~ SHPCRUISE + YEAR + MON + DAY
                      + LATITUDE + LONGITDE + T 
                      + GEAR + MESH + KRILL_BASIC + SST + CHLO + BATHY + MID_Z
                      + TIMEloc + day_of_year, 
                      data = krill, FUN = sum, na.action = NULL)

## For aggregated data where species are recorded as just "Euphausi" variant, we remove 
## data from the same sample that has more specific taxa information. We do this 
## assuming that these more speciated observations would be included in the basic
## "Euphausi" samples, and we don't want to record them twice
agg_krill <- agg_krill[order(agg_krill$KRILL_BASIC, decreasing = TRUE),] ## PUT ALL KRILL_BASIC OBSERVATIONS FIRST
non_krill_basic <- duplicated(agg_krill[,-c(13,19)]) ## IDENTIFY ALL NON-KRILL OBSERVATIONS FROM THE SAME
# SAMPLE AS COPEPODA AGGREGATE OBSERVATIONS
agg_krill <- agg_krill[c(!non_krill_basic),] ## REMOVE ALL FALSE (NON-COPEPODA REPLICATES)

krill <- agg_krill



library(stringr)
krill$UniqueSampleID <- paste0(krill$SHPCRUISE, krill$RECORD.ID) # Create a Unique ID for every net deployment
krill$UniqueSampleID <- str_sub(krill$UniqueSampleID, end=-6) # Remove the last 5 characters which are the 'taxaID'
krill$UniqueSampleID <- paste0(krill$UniqueSampleID, sprintf('%03.0f',krill$MID_Z)) # Now add in the depth as well

krill <- group_by(krill,UniqueSampleID)

krill_uni <- summarize(krill, 
                 Abundance = sum(VALUE.per.volu,na.rm = TRUE),
                 n_combined = n())
                 

krill <- filter(row_number() == 1)
                 
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
)


## PUT NA'S BACK
krill[krill$SST == -50, "SST"] <- NA
krill[krill$CHLO == -50, "CHLO"] <- NA
krill[krill$TIMEloc == -50, "TIMEloc"] <- NA

rm(agg_krill, non_krill_basic)

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
