library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)

imos <- read_rds("DatabaseOutput/IMOS_Final.rds")
sahfos <- read_rds("DatabaseOutput/SAHFOS_Final.rds")
copepod <- read_rds("DatabaseOutput/COPEPOD_Final.rds")
jedi <- read_rds("DatabaseOutput/JeDI_Final.rds")
L4 <- read_rds("DatabaseOutput/L4_Final.rds")
calcofi <- read_rds("DatabaseOutput/CalCOFI_Final.rds")
globec <- read_rds("DatabaseOutput/GLOBEC_Final.rds")
henschke <- read_rds("DatabaseOutput/Henschke_Final.rds")
warreen <- read_rds("DatabaseOutput/Warreen_Final.rds")
habas <- read_rds("DatabaseOutput/Habasque_Final.rds")
papa <- read_rds("DatabaseOutput/OSPapa_Final.rds")

dat <- rbind(copepod, jedi, imos, sahfos, L4, calcofi, globec, henschke, warreen, habas, papa)
# 
# dat <- rbind(globec, henschke, warreen)

dat2 <- dat %>% 
  mutate(
    Gear = as.numeric(Gear),
    Gear = replace(Gear, Project == "IMOS-NRS",1000),  # Code for NRS dropnet
    Gear = replace(Gear, Type == "CPR", 191), # Code for CPR
    Gear = replace(Gear, DBase == "JeDI", 299), # Gear code for JeDI. 299 is the first available gear number in COPEPOD
    Gear = replace(Gear, Project=="INSTOP-6", 300), # INSTOP-6 has missing gear - Create new gear Type
    Mesh = replace(Mesh, Type == "CPR", 270),
    Mesh = replace(Mesh, Project == "IMOS-NRS", 100),
    Upper_Z = replace(Upper_Z, Project == "Purcell_J", 1),
    Type = case_when(
      Gear == 191 ~ "CPR",
      Gear != 191 ~ "Net",
      is.na(Gear) ~ "Net"),
    Type = as.factor(Type),
    Gear_Mesh = str_c(as.character(Gear), as.character(Mesh), sep =  "_"), # Create Gear_Mesh re Heneghan et al
    Tow = case_when( # V = Vertical, O = Oblique, H = Horizontal, S = Surface (includes CPR)
      is.na(Tow) & Upper_Z == 0 ~ "V", # This is  a possibly dangerous assumption
      is.na(Tow) & Upper_Z != 0 ~ "H",
      Gear == 191 ~ "S",
      Project == "IMOS-NRS" ~ "V",
      !is.na(Tow) ~ Tow),
    Month = as.numeric(Month), # Create Month2 to standardise SH months to equivalent NH months
    Month2 = case_when(
      Latitude < 0 & Month <= 6 ~ Month + 6,
      Latitude < 0 & Month > 6 ~ Month - 6,
      Latitude >= 0 ~ Month),
    Date = ymd(sprintf('%04d%02d%02d', Year, Month, Day, sep = "")),
    DOY = yday(Date),
    DOY2 = DOY*0 - 999, # Calculate DOY2 (Day of Year, with SH coded as NH equivalent - i.e. summer SH coded same as summer NH)
    DOY2 = case_when(
      Latitude < 0 ~ DOY - 182,
      !is.na(Latitude) ~ DOY),
    DOY2 = case_when(
      DOY2 <= 0 ~ DOY2 + 365,
      DOY2 > 0 ~ DOY2),
    Upper_Z = case_when(
      Type == "CPR" ~ 0, # Include CPR data as 0-10 m 
      Project == "IMOS-NRS" ~ 0,
      !is.na(Upper_Z) ~ Upper_Z), # Include NRS data as 0-70 m
    Lower_Z = case_when(
      Type == "CPR" ~ 10,  # Include CPR data as 0-10 m 
      Project == "IMOS-NRS" ~ 70,
      !is.na(Lower_Z) ~ Lower_Z),
    Diff = Upper_Z - Lower_Z,
    Mid_Z = (Upper_Z + Lower_Z)/2,
    Gear = as.factor(Gear),
    Tow = as.factor(Tow),
    Gear_Mesh = as.factor(Gear_Mesh),
    ShipCruise = as.factor(ShipCruise),
    Transect = as.factor(Transect)) %>% 
  drop_na(TimeLocal) %>%  # Remove NAs
  drop_na(TotAbundance) %>% 
  filter(Diff <= 0) # Remove samples where Lower_Z deeper than Upper_Z

dat2 <- droplevels(dat2) 


# Add in Longhurst Provinces
Longhurst <- readOGR(dsn = "Longhurst", layer = "Longhurst_world_v4_2010")
Locations <- data.frame("Lat" = dat2$Latitude, "Lon" = dat2$Longitude, 
                        "Abund" = dat2$TotAbundance)
coordinates(Locations) <- ~Lon+Lat
crs(Locations) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
Provinces <- Longhurst[, "ProvCode"]
a.data <- over(Locations, Provinces)
dat2$Longhurst <- a.data$ProvCode
rm(a.data, Locations, Longhurst, Provinces)


# Check Abundance
summary(dat2$TotAbundance)
hist(log10(dat2$TotAbundance+1)) # 0 NA - Very high max abundance
dat2 <- dat2 %>% mutate(TotAbundance = replace(TotAbundance,TotAbundance>100000,100000))

# Check Gear
summary(dat2$Gear) # no MD
# Note that all JEDI has the same "Gear" (299)

# Check and code Type (CPR vs Net)
unique(dat2$Project[dat2$Type == "CPR"]) # Check - correct

# Check Mesh
hist(dat2$Mesh)
summary(dat2$Mesh) # No MD # Have confirmed 8000 is valid mesh size

summary(dat2$Tow) # No MD

# Uncorrected
hist(dat2$DOY)
summary(dat2$DOY)

hist(dat2$DOY2)
summary(dat2$DOY2)

# Check NewLocalTime
summary(dat2$TimeLocal)
hist(dat2$TimeLocal) # Very even

# Check Latitude and Longitude
summary(dat2$Latitude) # No MD
summary(dat2$Longitude) # No MD

hist(log10(abs(dat2$Diff)+1))

hist(log10(dat2$Lower_Z))

hist(dat2$Upper_Z)
summary(dat2$Upper_Z)

hist(dat2$Mid_Z)
summary(dat2$Mid_Z)

# Check Groups
table(dat2$Group)

# Check ShipCruise
levels(dat2$ShipCruise)

saveRDS(dat2, file = "DatabaseOutput/LatestDatabaseOuput_Final.rds")
######################## 4. END Pre-processing data #################
