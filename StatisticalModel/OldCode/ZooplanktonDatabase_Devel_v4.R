# ZooplanktonModels.R
# Ryan, Ant and Jase
# Last updated: 8/4/2019
# Feedback for Todd: 
# 1. Copepods - #/m2 converted to #/m3 by dividing by Water Strained, and it should be divide by net depth
# 2. Krill data has "-----" as UNITS
# 3. No zeros in the CPR data...

# Consider a Copepod only model (joining, carn and omni)


library(mgcv)
library(effects)
library(splines)
library(visreg)
library(stringr)
library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)

######## 0. COPEPOD data: Create Omnivorous and Carnivorous copepod groups ###############
# Only run Section 0 to recalculate Omnivores and Carnivores. Only run once to get .csv files
# 8 % of the records have unefined diets ebcause they are copepoda, or calenoida
#
# dat0 <- read.csv("copepods_2019.csv")
# # dat0 <- dat0 %>% rename(SpeciesName = SCIENTIFIC.NAME....modifiers...)
# 
# # Remove rows with NAs and null for VALUE.per.volu
# dat0 <- dat0 %>% filter(!is.na(VALUE.per.volu))
# 
# # Remove everything after -[ in SpeciesName
# Species_Namess <- sapply(strsplit(as.character(dat0$SCIENTIFIC.NAME....modifiers...), split = " -[", fixed = TRUE), function(x) (x[1]))
# dat0$SCIENTIFIC.NAME....modifiers... <- Species_Namess
# dat0$SCIENTIFIC.NAME....modifiers... <- as.factor(dat0$SCIENTIFIC.NAME....modifiers...)
# 
# # Import in FeedingType and clean it
# FeedingType <- read.csv("Copepod_Feeding_Type2.csv")
# 
# # Remove whitespace from end of FeedingType names
# FeedingType[, c(1, 2)] <- apply(FeedingType[, c(1, 2)], 2, 
#                                function(x){as.factor(str_trim(as.character(x), side = "right"))})
# FeedingType[, 1] <- as.factor(FeedingType[,1])
# 
# # Check species names the same in dat0$SpeciesName and FeedingType$TAXON_NAME - Yes!!!
# dplyr::setdiff(levels(FeedingType$TAXON_NAME), levels(dat0$SCIENTIFIC.NAME....modifiers...))
# 
# dat0 <- merge(dat0, FeedingType[, c("TAXON_NAME", "FEED")], 
#       by.x = "SCIENTIFIC.NAME....modifiers...", by.y = "TAXON_NAME")
# 
# # Set Diets: All copepod species assigned a diet
# levels(dat0$FEED)
# 
# # Diet: Carnivores (CC) vs Omnivores (CBF + CD + CH + CO + CP + SF)
# dat0$Diet <- "Omnivore"
# dat0$Diet[dat0$FEED == "CC"] <- "Carnivore"
# dat0$Diet <- as.factor(dat0$Diet)
# 
# # Diet2: Herbivore (CH) vs Carnivore (CC + CBF + CD + CH + CO + CP + SF)
# dat0$Diet2 <- "Carnivore"
# dat0$Diet2[dat0$FEED == "CH"] <- "Herbivore"
# dat0$Diet2 <- as.factor(dat0$Diet2)
# 
# # Change SpeciesName name back to original so same as in other files
# # dat0 <- dat0 %>% rename(SCIENTIFIC.NAME....modifiers... = SpeciesName)
# 
# datC <- dat0 %>% filter(Diet == "Carnivore")
# datC <- datC[, 1:(ncol(dat0) - 3)] # Remove last columns of diet information so same as other Group data files
# write.csv(datC, "copepods_carn_2019.csv", row.names = FALSE)
# 
# datO <- dat0 %>% filter(Diet == "Omnivore")
# datO <- datO[, 1:(ncol(dat0) - 3)] # # Remove last columns of diet information so same as other Group data files
# write.csv(datO, "copepods_omni_2019.csv", row.names = FALSE)  
# 
# rm(dat0, datC, datO, FeedingType, Species_Namess)

################ 1. COPEPOD data ################ 
################ 1A. COPEPOD data: Read in and select variables ################ 
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
dat6 <- read.csv("copepods_omni_2019.csv")
dat6 <- dat6 %>% mutate(Group = "OmniCopepods")
dat7 <- read.csv("copepods_carn_2019.csv")
dat7 <- dat7 %>% mutate(Group = "CarnCopepods")

dat <- rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7)
rm(dat1, dat2, dat3, dat4, dat5, dat6, dat7)

# Remove rows with NAs and null for VALUE.per.volu
dat <- dat %>% filter(!is.na(VALUE.per.volu))
summary(dat$VALUE.per.volu)
hist(log10(dat$VALUE.per.volu+1))

Retain <- c("Group", "SHPCRUISE", "YEAR", "MON", "DAY", "TIMEgmt", 
            "TIMEloc", "LATITUDE", "LONGITDE", "UPPER_Z", "LOWER_Z", 
            "T", "GEAR", "MESH", "Water.Strained", "Original.VALUE", "Orig.UNITS", "VALUE.per.volu", "UNITS", 
            "SCIENTIFIC.NAME....modifiers...", "RECORD.ID", "DATASET.ID", "SHIP", "INST", "SST", "CHLO", "BATHY")
dat <- dat[, Retain]
dat <- dat %>% rename(Group = Group, ShipCruise = SHPCRUISE, Year = YEAR, Month = MON, 
                      Day = DAY, TimeGMT = TIMEgmt, TimeLoc = TIMEloc, 
                      Latitude = LATITUDE, Longitude = LONGITDE, Upper_Z = UPPER_Z, 
                      Lower_Z = LOWER_Z, Tow = T, Gear = GEAR, Mesh = MESH, 
                      Volume = Water.Strained, OrigValue = Original.VALUE, OrigUnits = Orig.UNITS,
                      Abundance = VALUE.per.volu, Units = UNITS, SpeciesName = SCIENTIFIC.NAME....modifiers..., 
                      Record.ID = RECORD.ID, Project = DATASET.ID, Ship = SHIP, 
                      Inst = INST, SST = SST, Chl = CHLO, Bathy = BATHY)

################ 1B. COPEPOD data: Remove CPR double counting ################ 
# Remove all inverted commas within Species Name - causing problems with str_detect below
dat <- dat %>% mutate_at("SpeciesName", str_replace_all, 
                         pattern = '["]', replacement = "")
dat$SpeciesName <- as.factor(dat$SpeciesName)

# Split the dataset into two so that we can check the CPR data
# Ant and Jase could not find a good dplyr way to do this in one step
CPR <- filter(dat, Project == "SAHFOS-CPR Atlantic Ocean")
Other <- filter(dat, Project != "SAHFOS-CPR Atlantic Ocean")

# filter out the double-counting in the CPR dataset. The \\ cancels the special characters for regexp
# Just remove duplicates. Keep all eggs and nauplii from all Groups - will be low abundance
CPR <- CPR %>% dplyr::filter(!str_detect(SpeciesName, "Calanus CPR\\-TOTAL traverse") &
                               !str_detect(SpeciesName, "Calanus TOTAL \\(total count\\)") &
                               !str_detect(SpeciesName, "Centropages chierchiae CPR\\-TOTAL traverse") &
                               # !str_detect(SpeciesName, "Copepoda \\-\\[ eggs \\]\\-") &
                               !str_detect(SpeciesName, "Copepoda TOTAL \\(total count\\) \\-\\[ \\]\\-") &
                               !str_detect(SpeciesName, "Metridia CPR\\-TOTAL traverse") &
                               !str_detect(SpeciesName, "Euphausiacea TOTAL") &
                               # !str_detect(SpeciesName, "Euphausiacea \\-\\[ eggs \\]\\-") &
                               !str_detect(SpeciesName, "Chaetognatha CPR\\-TOTAL traverse") &
                               !str_detect(SpeciesName, "Thaliacea"))

dat <- bind_rows(CPR, Other)
rm(CPR, Other)

######## 1C. COPEPOD data: Calculate TotAbundance for each Group per Sample ###############
# Unique ID for each sample
dat <- dat %>% mutate(UniqueSampleID = paste0(ShipCruise, Record.ID), # Create Unique ID for every net deployment
                      UniqueSampleID = str_sub(UniqueSampleID, end = -6), # Remove last 5 characters which are the 'taxaID')
                      UniqueSampleID = paste0(UniqueSampleID, sprintf('%04.0f', Upper_Z)), # Now add in depth as well
                      UniqueSampleID = as.factor(UniqueSampleID)) 

# Arrange so by Group and UniqueSampleID
dat <- dat %>% arrange(Group, UniqueSampleID)

# Calculate total abundance for each UniqueSampleID
dat2 <- dat %>% mutate(Group = as.factor(Group),
                       UniqueSampleID = as.factor(UniqueSampleID)) %>%
                       group_by(Group, UniqueSampleID) %>%
                       summarise(TotAbundance = sum(Abundance))

# All rows in a UniqueSampleID are the same except for ScientificName, so keep the distinct ones 
d <- dat %>% mutate(Group = as.factor(Group)) %>%
  group_by(Group, UniqueSampleID) %>%
  distinct(Group, UniqueSampleID, .keep_all = TRUE)
dat <- right_join(d, dat2, by = c("Group", "UniqueSampleID"))
rm(d, dat2) # Clean up
dat <- dat %>% dplyr::select(-Abundance)

######## 1D. COPEPOD data: Add zeros to CPR data ###############
dat <- ungroup(dat) # Make sure the df is ungrouped
CPR <- filter(dat, Project == "SAHFOS-CPR Atlantic Ocean") # Split into CPR and other
Other <- filter(dat, Project != "SAHFOS-CPR Atlantic Ocean")

CPR_wide <- spread(CPR, Group, TotAbundance, fill = 0) # Wide format: Move the Abundance into each group. 
CPR_wide <- arrange(CPR_wide, ShipCruise, Year, Month, Day, TimeGMT) # Arrange to make it easier to debug
CPR_wide <- group_by(CPR_wide, UniqueSampleID) # Regroup by UniqueSampleID

# Calculate total abundance for each UniqueSampleID
CPR_wide2 <- summarise(CPR_wide,
                       Chaetognaths = sum(Chaetognaths),
                       CarnCopepods = sum(CarnCopepods),
                       OmniCopepods = sum(OmniCopepods),
                       Euphausiids = sum(Euphausiids),
                       Jellyfish = sum(Jellyfish),
                       Larvaceans = sum(Larvaceans),
                       Salps = sum(Salps))
# Remove the Groups from CPR_wide so we don't get duplicates on the inner_join
CPR_wide <- dplyr::select(CPR_wide,-c("Chaetognaths", "CarnCopepods", "OmniCopepods", "Euphausiids", "Jellyfish", "Larvaceans", "Salps"))
CPR_wide3 <- inner_join(CPR_wide, CPR_wide2, by="UniqueSampleID") # Join CPR_wide (metadata) with CPR_wide2 (Abbundance)

CPR_long <- gather(CPR_wide3, Group, TotAbundance, Chaetognaths:Salps, factor_key = TRUE)

CPR_long <- arrange(CPR_long, ShipCruise, Year, Month, Day, TimeGMT) # Arrange to make it easier to debug
CPR_long <- distinct(CPR_long, Group, UniqueSampleID, .keep_all = TRUE)
dat <- bind_rows(CPR_long, Other)
rm(CPR, CPR_long, CPR_wide, CPR_wide2, CPR_wide3, Other)
dat <- dat %>% dplyr::select(-SpeciesName)

#########################
## FINISH ADDING ZEROS ##
#########################

######## 1E. COPEPOD data: Final clean up ###############
# Only include Original units that are #/m3 and #/haul and drop dodgy ----- as Units
dat <- dat %>% filter(!str_detect(Units, "-----")) %>%
  filter(str_detect(OrigUnits, "#/m3") | str_detect(OrigUnits, "#/haul")) %>%
  droplevels() %>%
  mutate(OrigUnits = fct_recode(OrigUnits, "#/m3" = "      #/m3", "#/haul" = "    #/haul"), 
         Units = fct_recode(Units, "#/m3" = " #/m3"))

# Clean up time of day
hist(dat$TimeLoc)
hist(dat$TimeGMT)
dat <- dat %>% mutate(TimeLocal = case_when(
  TimeLoc >= 0 ~ TimeLoc,
  TimeGMT >= 0 & TimeGMT <= 24 & TimeLoc == -99 ~ TimeGMT + Longitude/15)) 
dat <- dat %>% mutate(TimeLocal = case_when(
  TimeLocal >= 0 & TimeLocal <= 24 ~ TimeLocal,
  TimeLocal > 24 ~ TimeLocal - 24,
  TimeLocal < 0 ~ TimeLocal + 24))
hist(dat$TimeLocal)
dat <- dat %>% dplyr::select(-c("TimeLoc"))

dat$DBase <- as.factor("COPEPOD") # Assign level to databbase origin
dat$ShipCruise <- as.factor(dat$ShipCruise)
dat$Project <- as.factor(dat$Project)
dat$Group <- as.factor(dat$Group)

################ 2. IMOS data ################ 
dat2 <- read.csv("imos_2019v3.csv")

# Recode Larvacean to Larvaceans, Chaetognath to Chaetognaths, Krill to Euphausiids
dat2 <- gather(dat2, Group, Abundance, Larvacean:Salps, factor_key = TRUE)
dat2 <- dat2 %>% mutate(Group = fct_recode(Group, "Larvaceans" = "Larvacean", 
                                           "Chaetognaths" = "Chaetognath", 
                                           "Euphausiids" = "Krill",
                                           "CarnCopepods" = "Carnivore.Cope",
                                           "OmniCopepods"= "Ominvore.Cope"))

dat2 <- dat2 %>% rename(Project = PROJECT, ShipCruise = SHPCRUISE, Latitude = LATITUDE, Longitude = LONGITDE, 
                        Year = YEAR, Month = MON, Day = DAY, TimeGMT = TimeGMT, SST = SST, Chl = CHLO, 
                        Bathy = BATHY, Group = Group, TotAbundance = Abundance)

imos_time <- hms(dat2$TimeGMT)
dat2$TimeGMT <- hour(imos_time) + (minute(imos_time) + (second(imos_time)/60))/60
rm(imos_time)

# Should we be rounding this Longitude/15 to the nearest integer to represent timezones?
dat2$TimeLocal <- dat2$TimeGMT + dat2$Longitude/15 # Add 1hr for every 15 degrees of longitude

dat2 <- dat2 %>% mutate(TimeLocal = case_when(
  TimeLocal >= 0 & TimeLocal <= 24 ~ TimeLocal,
  TimeLocal > 24 ~ TimeLocal - 24,
  TimeLocal < 0 ~ TimeLocal + 24))
hist(dat$TimeLocal)

dat2 <- dat2 %>% dplyr::select(-TimeGMT)

dat2 <- dat2 %>% mutate(Project = as.character(Project), 
  Type = case_when(
  Project == "cpr" | Project == "SO" ~ "CPR",
  Project == "nrs" ~ "Net"),
  Project = case_when(
    Project == "cpr" ~ "IMOS-CPR",
    Project == "nrs" ~ "IMOS-NRS",
    Project == "SO" ~ "SO-CPR"),
  Type = as.factor(Type), 
  Project = as.factor(Project),
  Group = as.factor(Group))

dat2$DBase <- as.factor("IMOS") # Assign level to databbase origin

# Bind data
dat <- bind_rows(dat, dat2)
rm(dat2)

################ 3. JEDI data ################ 
jedi <- read.csv("jedi_2019_FullDataset_Enviro.csv",na.strings = "NA")

jedi[jedi==-9999] <- NA # jedi missing data indicator

# jedi$day <- as.numeric(levels(jedi$day))[jedi$day] # Convert to numeric
# jedi$depth_upper <- as.numeric(levels(jedi$depth_upper))[jedi$depth_upper] # Convert to numeric
# jedi$depth_lower <- as.numeric(levels(jedi$depth_lower))[jedi$depth_lower] # Convert to numeric


# Retain <- c("taxon", "rank_phylum", "year" , "month", "day", "time_local", "lat",
#             "lon", "depth_upper", "depth_lower", "net_mesh", "collection_method",
#             "VALUE.per.volu","owner_dataset","project_title",
#             "SST", "CHLO", "BATHY","data_type")

Retain <- c("taxon", "rank_phylum", "year" , "month", "day", "time_local", "lat",
            "lon", "depth", "depth_upper", "depth_lower", "net_mesh", "collection_method",
            "density","owner_dataset","project_title",
            "SST", "CHLO", "BATHY","data_type")

jedi <- jedi %>% dplyr::select(Retain)
jedi <- jedi %>% rename(Year = year, Month = month,
                        Day = day, TimeLoc = time_local,
                        Latitude = lat, Longitude = lon, RecordedDepth = depth, Upper_Z = depth_upper,
                        Lower_Z = depth_lower, Gear = collection_method, Mesh = net_mesh,
                        Abundance = density,  Project = owner_dataset,
                        Inst = project_title,
                        SST = SST, Chl = CHLO, Bathy = BATHY, data_type=data_type)

jedi <- jedi[jedi$Gear == "plankton_net",] # Only quantitative surveys
jedi <- jedi[jedi$data_type == "quantitative",] # Only quantitative surveys
# jedi <- jedi[jedi$Mesh != 'nd',] # Only surveys with mesh size recorded
jedi <- jedi %>% filter(!is.na(Mesh)) # Filter for NA
jedi <- jedi[jedi$Inst != 'COPEPOD',] # No surveys from the COPEPOD database

# jedi$TotAbundance <- as.numeric(levels(jedi$TotAbundance))[jedi$TotAbundance] # Convert to numeric
jedi <- jedi %>% filter(!is.na(Abundance)) # Filter for NA

jedi$rank_phylum[jedi$taxon == "medusa"] <- "Cnidaria" # Generate some phyla to make things easier below
jedi$rank_phylum[jedi$taxon == "jellies"] <- "Cnidaria"

# Only keep the required phyla
jedi <- filter(jedi, taxon == "salp" | rank_phylum == 'Cnidaria')

# Move the relevent info into Group
jedi <- jedi%>% add_column(Group=NA)
jedi$Group[jedi$taxon=="salp"] <- "Salps"
jedi$Group[jedi$rank_phylum=="Cnidaria"] <- "Jellyfish"
jedi$Group <- as.factor(jedi$Group)

# Unique ID for each sample
jedi <- jedi %>% mutate(UniqueSampleID = paste0(sprintf('%02.2f', abs(Latitude)),
                                                sprintf('%03.2f', abs(Longitude)),
                                                sprintf('%4d', Year),
                                                sprintf('%2d', Month),
                                                sprintf('%2d', Day)), # Create Unique ID for every net deployment
                      UniqueSampleID = as.factor(UniqueSampleID)) # Now add in depth as well

# Arrange so by Group and UniqueSampleID
jedi <- jedi %>% arrange(Group, UniqueSampleID)

# Calculate total abundance for each UniqueSampleID
jedi2 <- jedi %>% mutate(Group = as.factor(Group),
                         UniqueSampleID = as.factor(UniqueSampleID)) %>%
                           group_by(Group, UniqueSampleID) %>%
                           summarise(TotAbundance = sum(Abundance))
jedi <- ungroup(jedi)

# All rows in a UniqueSampleID are the same except for ScientificName, so keep the distinct ones
d <- jedi %>% mutate(Group = as.factor(Group)) %>%
  group_by(Group, UniqueSampleID) %>%
  distinct(Group, UniqueSampleID, .keep_all = TRUE)
jedi <- right_join(d, jedi2, by = c("Group", "UniqueSampleID"))

rm(d, jedi2) # Clean up

# jedi <- dplyr::select(jedi,-Abundance)


# Add missing variables
jedi <- jedi%>% add_column(Volume=NA, OrigValue=NA, OrigUnits="m3", Units="m3",
                           Tow=NA, ShipCruise=NA, TimeGMT=NA, Type="Net", Record.ID=NA, Ship=NA)

jedi$ShipCruise=jedi$Inst # Use project as the Ship-Cruise random variable
jedi$UniqueSampleID=jedi$Inst # Use project as the UniqueSampleID random variable
jedi <- droplevels(jedi)

# Mesh is in mm. Convert to um
jedi$Mesh <- as.numeric(jedi$Mesh)*1e3

# Some incorrect meshes in the database
jedi$Mesh[jedi$Inst=="Helgoland_Roads"] <- 500 # https://www.researchgate.net/publication/226735944_Helgoland_Roads_North_Sea_45_Years_of_Change
jedi$Mesh[jedi$Inst=="PACES"] <- 76 # https://www.researchgate.net/publication/226910458_Martens_P_van_Beusekom_JEE_Zooplankton_response_to_a_warmer_northern_Wadden_Sea_Helgol_Mar_Res_62_67-75

# Sample from the surface. Adding in Upper/Lower
jedi$Upper_Z[jedi$Inst=="PACES"] <- 0
jedi$Lower_Z[jedi$Inst=="PACES"] <- 0

# Some time are numeric
jedi_time <- hms(jedi$TimeLoc) # Some fail due to na
jedi$TimeLocal <- hour(jedi_time) + (minute(jedi_time) + (second(jedi_time)/60))/60
rm(jedi_time)
hist(jedi$TimeLocal)

# Remove extra variables
jedi <- jedi %>% dplyr::select(-c(taxon, rank_phylum, data_type, TimeLoc, taxon, rank_phylum))

jedi$DBase <- as.factor("JeDI") # Assign level to databbase origin

dat <- bind_rows(dat,jedi)
dat <- droplevels(dat) # Make sure all the extra levels are gone.

saveRDS(dat, file = "LatestDatabaseOuput_NoQC.rds")
######################## 4. Pre-processing data #################

dat <- readRDS("LatestDatabaseOuput_NoQC.rds")

# Check Abundance
summary(dat$TotAbundance)
hist(log10(dat$TotAbundance+1)) # 1 NA - looks OK
dat <- dat %>% filter(!is.na(TotAbundance)) # Remove NAs

# Check Gear
summary(dat$Gear)
dat$Gear[dat$Project == "nrs"] <- 1000 # Code for NRS dropnet
dat$Gear[dat$Type == "CPR"] <- 191 # Code for CPR
dat$Gear[dat$DBase == "JeDI"] <- 299 # Gear code for JeDI. 299 is the first available gear number in COPEPOD
dat$Gear[dat$Project=="INSTOP-6"] <- 300 # INSTOP-6 has missing gear - I will create a new gear type
dat$Gear[dat$DBase=="IMOS" & dat$Type=="Net"] <- 301 # Add in gear for NRS
summary(dat$Gear) # no MD
# Note that all JEDI has the same "Gear" (299)

# Check and code Type (CPR vs Net)
dat <- dat %>% mutate(Type = case_when(
  Gear == 191 ~ "CPR",
  Gear != 191 ~ "Net",
  is.na(Gear) ~ "Net"))
dat$Type <- as.factor(dat$Type)
unique(dat$Project[dat$Type == "CPR"]) # Check - correct

# Check Mesh
hist(dat$Mesh)
summary(dat$Mesh)
sum(dat$Mesh == -999.0)
dat$Mesh[dat$Type == "CPR"] <- 270
dat$Mesh[dat$Project == "IMOS-NRS"] <- 100
dat <- dat %>% filter(Mesh != -999.0) # Remove rows with no mesh size given
summary(dat$Mesh) # No MD # Have confirmed 8000 is valid mesh size

## Create Gear_Mesh variable: Heneghan et al
dat$Gear_Mesh <- as.factor(paste(as.character(dat$Gear), as.character(dat$Mesh), sep =  "_"))
levels(dat$Gear_Mesh)

# Check Tow
# V = Vertical, O = Oblique, H = Horizontal, S = Surface (includes CPR)
summary(dat$Tow)
dat$Tow[dat$Project == "IMOS-NRS"] <- "V" # Code for NRS dropnet
dat$Tow[dat$Type == "CPR"] <- "S" # Code for CPR
dat$Tow[dat$Tow == "-" & dat$Upper_Z == 0] <- "V"
dat$Tow[dat$Tow == "-" & dat$Upper_Z != 0] <- "H"
dat$Tow[is.na(dat$Tow) & dat$Upper_Z == 0] <- "V"
dat$Tow[is.na(dat$Tow) & dat$Upper_Z != 0] <- "H"
dat <- droplevels(dat)
summary(dat$Tow)
# upper_z <- dat[is.na(dat$Upper_Z),]


# Check and calculate Months and DOY
# Create  variable (Month2) to standardise SH months to equivalent NH months
dat <- dat %>% mutate(Month = as.numeric(Month), 
                      Month2 = case_when(
                        Latitude < 0 & Month <= 6 ~ Month + 6,
                        Latitude < 0 & Month > 6 ~ Month - 6,
                        Latitude >= 0 ~ Month))

# Calculate DOY (Day of Year)
dat$Date <- ymd(sprintf('%04d%02d%02d', dat$Year, dat$Month, dat$Day, sep = ""), tz = "Europe/London")
dat$DOY <- yday(dat$Date)
hist(dat$DOY)
summary(dat$DOY)

# Calculate DOY2 (Day of Year, with SH coded as NH equivalent - i.e. summer SH coded same as summer NH)
dat$DOY2 <- dat$DOY*0 - 999
dat$DOY2[dat$Latitude<0] <- dat$DOY[dat$Latitude<0] - 182
dat$DOY2[dat$DOY2<0] <- dat$DOY2[dat$DOY2<0] + 365
dat$DOY2[dat$Latitude>=0] <- dat$DOY[dat$Latitude>=0]
hist(dat$DOY2)
summary(dat$DOY2)

# Check NewLocalTime
summary(dat$TimeLocal)
dat <- dat %>% filter(!is.na(TimeLocal)) # Remove NAs
hist(dat$TimeLocal)

# Check Latitude and Longitude
summary(dat$Latitude) # No MD
summary(dat$Longitude) # No MD

# Check Bathy - Note: 0s must be where it is shallow - recode them to -20 m
# Convert all bathymetry to positive
summary(dat$Bathy) # No missing data
hist(dat$Bathy)
sum(dat$Bathy < 0)
sum(dat$Bathy == 0)
sum(dat$Bathy > 0)
dat$Bathy[dat$Bathy == 0] <- -20
dat <- dat %>% mutate(Bathy = abs(Bathy), # Make Bathy positive
                      Bathy = replace(Bathy, Bathy > 7000, 7000)) # Set max bathymetry
hist(dat$Bathy)


# Check Lower_Z and Upper_Z
# Include CPR data as 0-10 m and NRS data as 0-60 m
dat$Upper_Z[dat$Type == "CPR"] <- 0
dat$Lower_Z[dat$Type == "CPR"] <- 10

dat$Upper_Z[dat$Project == "IMOS-NRS"] <- 0
dat$Lower_Z[dat$Project == "IMOS-NRS"] <- 70

dat <- dat %>% mutate(Diff = Upper_Z - Lower_Z)
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
summary(dat$Mid_Z)

# Check SST
summary(dat$SST) # 7807 NAs
dat <- dat %>% filter(!is.na(SST)) # Remove NAs
hist(dat$SST)

# Check Chl
summary(dat$Chl) # 10467 NAs
dat <- dat %>% filter(!is.na(Chl)) # Remove NAs
hist(dat$Chl)
hist(log10(dat$Chl))
dat <- dat %>% mutate(Chl = replace(Chl, Chl > 10, 10)) # Set max Chl

# Check Groups
dat$Group <- as.factor(dat$Group)
table(dat$Group)

# Check ShipCruise
dat$ShipCruise <- as.factor(dat$ShipCruise)
levels(dat$ShipCruise)

# Add in Longhurst Provinces
Longhurst <- readOGR(dsn = "Longhurst", layer = "Longhurst_world_v4_2010")
# curr_zoo <- krill
Locations <- data.frame("Lat" = dat$Latitude, "Lon" = dat$Longitude, 
                        "Abund" = dat$TotAbundance)
coordinates(Locations) <- ~Lon+Lat
crs(Locations) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
Provinces <- Longhurst[, "ProvCode"]
a.data <- over(Locations, Provinces)
dat$Longhurst <- a.data$ProvCode
rm(a.data, Locations, Longhurst, Provinces)

# Latitude cut-off - e.g. 15% of data beyond this Latitude
hist(dat$Latitude)
Percent <- 2.5 # In each tail
Length <- length(dat$Latitude)
# Length <- 100
LatSorted <- sort(dat$Latitude)
Lower <- round(Percent / 100 * Length)
Upper <- round((100 - Percent) / 100 * Length)
LatLower <- LatSorted[Lower]
LatUpper <- LatSorted[Upper]
# Run this line to use a lat cut-off
# dat <- dat %>% filter(Latitude >= LatLower & Latitude <= LatUpper)
rm(LatSorted)

dat$Gear <- as.factor(dat$Gear)
dat$Project <- as.factor(dat$Project)
dat$DBase <- as.factor(dat$DBase)

saveRDS(dat, file = "LatestDatabaseOuput_Final.rds")
######################## 4. END Pre-processing data #################
