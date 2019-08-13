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

work_direct_main ="~/Library/Mobile Documents/com~apple~CloudDocs/PhD/" # Set your main directory where you
# put "Zooplankton GAMs"

work_direct = paste(work_direct_main, "Zooplankton GAMs/Larvaceans", sep = "")

setwd(work_direct)

larvaceans <- within(read.csv("larvs_2019.csv"),{
  
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

## PLOT THE DATA - THIS TAKES A FEW MINUTES
#source(paste(work_direct_main, "Zooplankton GAMs/Zoo Data Summary Plots.R", sep = ""))
#plot_zoo(larvaceans, "larvs_2019", save_plots = FALSE)


## REMOVE CPR DATA
#cpr_data <- grep("CPR", larvaceans$DATASET.ID)
#if(sum(cpr_data) > 0){
#larvaceans <- larvaceans[-cpr_data,]
#}

## REMOVE ALL ORIG.VALUE THAT ISN'T # M^-3
larvaceans <- larvaceans[which(larvaceans$Orig.UNITS %in% c("      #/m3", "    #/haul")),]

## Remove unecessary columns
larvaceans <- larvaceans[,-c(14:21,23,26:35,37:42)] # unnecessary columns

# Remove presence/absence points
larvaceans <- larvaceans[!is.na(larvaceans$TOT_ABUND),] #  presence/absence points

# REMOVE MESH SIZES OUTSIDE 100-500UM ()
#larvaceans <- larvaceans[c(larvaceans$MESH >= 00 & larvaceans$MESH <= 500),] #  mesh above 500um (~30 obs)
###

## REMOVE SAMPLES FROM BEFORE 1958
larvaceans <- larvaceans[c(larvaceans$YEAR >= 1958),]

## CLEAN UP SPECIES NAMES
species_names <- as.character(larvaceans$SCIENTIFIC.NAME....modifiers...)
species_namess <- sapply(strsplit(species_names, split = " -[", fixed = TRUE), function(x) (x[1]))
larvaceans$SCIENTIFIC.NAME....modifiers...<- species_namess
colnames(larvaceans)[17] <- "GROUP"

################# IMOS DATA ##########################
work_direct_imos = paste(work_direct_main, "Zooplankton GAMs/imos_2019.csv", sep = "") # OR "Zooplankton GAMs/imos_2.csv if you want to use 2017 imos data
imos <- read.csv(work_direct_imos)
curr_zoo <- "Larvacean"
wanted <- c("SHPCRUISE", "LATITUDE", "LONGITDE", "YEAR", "MON", "DAY", "TIMEloc", curr_zoo, "BATHY", "CHLO", "SST")
imos <- imos[,c(wanted)]
colnames(imos)[8] <- "TOT_ABUND"
imos$MESH <- 100
imos$GEAR <- 191 # CODED AS 'CPR' GEAR FOR ALL IMOS OBS - COPEPOD DATABASE GEAR NUMBER 191
imos$GROUP <- as.character("Larvacean")
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

larvaceans$SHPCRUISE <- as.character(larvaceans$SHPCRUISE)

larvaceans <- merge(larvaceans, imos, all = TRUE)

larvaceans$SHPCRUISE <- as.factor(larvaceans$SHPCRUISE)

########## IMPORT JEDI DATA
work_direct_jedi = paste(work_direct_main, "Zooplankton GAMs/jedi_2019.csv", sep = "") # OR "Zooplankton GAMs/imos_2.csv if you want to use 2017 imos data
jedi <- read.csv(work_direct_jedi)

tt = jedi[jedi$taxon == 'larvacean',] # Only larvaceans
tt = tt[tt$data_type == 'quantitative',] # Only quantitative surveys
tt = tt[tt$net_mesh != 'nd',] # Only surveys with mesh size recorded
tt = tt[tt$project_title != 'COPEPOD',] # No surveys from the COPEPOD database

jedi = tt

# Mean sample depth "MID_Z"
jedi$MID_Z <- (as.numeric(as.character(jedi$depth_upper)) + as.numeric(as.character(jedi$depth_lower)))/2 

## Rename variables to match COPEPOD file and pull out variables we want
jedi <- jedi %>% rename(SHPCRUISE = project_title, YEAR = year, DAY = day, TIMEloc = time_local,
                        TOT_ABUND = VALUE.per.volu, MESH = net_mesh, GEAR = collection_method, GROUP = taxon)
wanted <- c("SHPCRUISE", "LATITUDE", "LONGITDE", "YEAR", "MON", "DAY", "TIMEloc", "TOT_ABUND", 
            "BATHY", "CHLO", "SST", "MESH", "GEAR", "GROUP", "MID_Z")
jedi<- jedi[,c(wanted)]

jedi$T <- NA # Tow type not recorded
jedi$MON <- as.numeric(as.character(jedi$MON))
jedi$day_of_year <- round((as.integer((jedi$MON))-1)*30.4 + as.integer(as.character(jedi$DAY)))

# Standardise Southern Hemisphere months to equivalent Northern Hemisphere months
jedi$SH_1 <- c(jedi$LATITUDE < 0 & jedi$MON <=6) # Jan - June in South Hempishphere
jedi$SH_2 <- c(jedi$LATITUDE < 0 & jedi$MON  > 6) # July - December in South Hemisphere
jedi$MON[jedi$SH_1] <- jedi$MON[jedi$SH_1] + 6
jedi$MON[jedi$SH_2] <- jedi$MON[jedi$SH_2] - 6

jedi$BATHY <- abs(jedi$BATHY)

larvaceans$SHPCRUISE <- as.character(larvaceans$SHPCRUISE)

larvaceans <- merge(larvaceans, jedi, all = TRUE)

####### AGGREGATE speciated SAMPLES
larvaceans[is.na(larvaceans$CHLO),"CHLO"] <- -50
larvaceans[is.na(larvaceans$SST),"SST"] <- -50
larvaceans[is.na(larvaceans$TIMEloc), "TIMEloc"] <- -50

agg_larvaceans <- aggregate(TOT_ABUND ~ SHPCRUISE + YEAR + MON + DAY + LATITUDE + 
                              LONGITDE + T + GEAR + MESH + SST + CHLO + BATHY +
                              MID_Z + TIMEloc + day_of_year, 
                         data = larvaceans, FUN = sum, na.action = NULL)

agg_larvs <- agg_larvs[order(agg_larvs$APPENDICULARIA, decreasing = TRUE),] ## PUT ALL APPENDICULARIA OBSERVATIONS FIRST

non_appendicularia <- duplicated(agg_larvs[,-c(13,19)]) ## IDENTIFY ALL APPENDICULARIA OBSERVATIONS FROM THE SAME
## SAMPLE AS APPENDICULARIA OBSERVATIONS
agg_larvs <- agg_larvs[c(!non_appendicularia),] ## REMOVE ALL FALSE (NON-APPENDICULARIA REPLICATES)

larvaceans <- agg_larvaceans

## PUT NA'S BACK
larvaceans[larvaceans$SST == -50, "SST"] <- NA
larvaceans[larvaceans$CHLO == -50, "CHLO"] <- NA
larvaceans[larvaceans$TIMEloc == -50, "TIMEloc"] <- NA

larvaceans$GEAR_MESH <- as.factor(paste(as.character(larvaceans$GEAR), as.character(larvaceans$MESH), sep = "_"))
larvaceans$log10CHLO <- log10(larvaceans$CHLO)

## ADD IN LONGHURST PROVINCES
setwd(paste(work_direct_main, "Zooplankton GAMs", sep = ""))
longhurst <- readOGR(dsn = "Longhurst", layer = "Longhurst_world_v4_2010")
curr_zoo <- larvaceans
curr_locs <- data.frame("lat" = curr_zoo$LATITUDE, "lon" = curr_zoo$LONGITDE, 
                        "abund" = curr_zoo$TOT_ABUND)
coordinates(curr_locs) <- ~lon+lat
crs(curr_locs) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
provinces <- longhurst[,"ProvCode"]
a.data <- over(curr_locs, provinces)
curr_zoo$LONGHURST <- a.data$ProvCode
setwd(work_direct)

larvaceans <- curr_zoo

#### CLEAN UP TIMES, SO THAT THEY ARE ALL NUMERIC
curr_zoo <- larvaceans

fixtimes <- grep(":", curr_zoo$TIMEloc) # Which times are recorded as "HH:MM:SS"
fixtimes2 <- curr_zoo[fixtimes, "TIMEloc"] # Pull out these times

fixedtimes <- sapply(strsplit(fixtimes2,":"), # Convert to numeric equivalent
                     function(x) {
                       x <- as.numeric(x)
                       x[1]+x[2]/60
                     }
)

curr_zoo[fixtimes, 'TIMEloc'] <- as.character(fixedtimes)
larvaceans <- curr_zoo

write.csv(larvaceans, "clean_larvs.csv", row.names = FALSE)

## PLOT THE DATA - THIS TAKES A FEW MINUTES
#source(paste(work_direct_main, "Zooplankton GAMs/Zoo Data Summary Plots.R", sep = ""))
#plot_zoo(larvaceans, "larvaceans")

#####################################################################
#####################################################################
larvaceans <- read.csv("clean_larvs.csv")
larv_t <- larvaceans#[!is.na(larvaceans$TIMEloc),] ## REMOVE NON TIME-STAMP OBSERVATIONS
larv_t[larv_t$MESH < 1000, 'MESH'] <- 1000 ## CONSTRAIN MESH SIZES TO <= 1000UM
larv_t <- larv_t[larv_t$MID_Z <= 200,] ## REMOVE OBS FROM OVER 200M DEPTH

############################################################
## LAT CUT OFF - WHERE 15% OF THE DATA IS AFTER THIS LATITUDE
lat_cut = data.frame("lat" = seq(0, 90, 1), "num" = NA)
curr_zoo = larv_t
curr_zoo$latits = abs(curr_zoo$LATITUDE)

for(i in 1:dim(lat_cut)[1]){
  curr_lat = lat_cut[i,1]
  lat_cut[i,"num"] =  sum((curr_zoo$latits >= curr_lat), na.rm = TRUE) 
}

lat_cut = lat_cut[min(which(lat_cut[1:90,"num"] < 0.15*dim(curr_zoo)[1])), "lat"]

curr_zoo = curr_zoo[c(abs(curr_zoo$LATITUDE) < lat_cut),]
###################################################

larv_t = curr_zoo ## If you want to keep all the latitudes, don't run this line

min_val = min(larv_t[larv_t$TOT_ABUND > 0, 'TOT_ABUND'], na.rm =TRUE)/2

## LARVACEANS MODEL 1
gm1 <- gam(log10(TOT_ABUND+min_val) ~   s(log10CHLO, fx = T, k = 5) +
             te(SST, day_of_year, bs = c("cr", "cc")) + s(BATHY, fx = T, k = 5) 
           + GEAR_MESH + s(TIMEloc), 
           data = larv_t, na.action = na.exclude)
summary(gm1)

# PLOT RESIDUALS
par(mfrow = c(2,2))
gam.check(gm1)
abline(0,1, lwd = 2, col = "red")
par(mfrow = c(1,1))

## SAVE GAM
gam_name = "larvs_gam11.rds"
gam_save_file = paste(work_direct_main, "Zooplankton GAMs/Interpolate GAMs/", gam_name, sep = "")
saveRDS(gm1, file = gam_name) 
saveRDS(gm1, file = gam_save_file) # Save to Interpolate GAMs file

## SAVE PLOTS
png(filename="larv_gam.png", units="in",  width=7, height=6, res=600)
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
