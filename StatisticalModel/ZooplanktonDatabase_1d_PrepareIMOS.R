rm(list = ls())

library(tidyverse)
library(lubridate)
library(stringr)

Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")

################ 2. IMOS data ################ 
imos <- read.csv("Data/IMOSv4.csv")

# Recode Larvacean to Larvaceans, Chaetognath to Chaetognaths, Krill to Euphausiids
imos <- imos %>% rename(Larvaceans = Larvacean, 
                        Chaetognaths = Chaetognath, 
                        Euphausiids = Krill,
                        CarnCopepods = Carnivore.Cope,
                        OmniCopepods = Ominvore.Cope,
                        Year = YEAR, Month = MON, Day = DAY,
                        Project = PROJECT, ShipCruise = SHPCRUISE) %>% 
  gather(c("Chaetognaths", "CarnCopepods", "OmniCopepods", "Euphausiids", "Jellyfish", "Larvaceans", "Salps"), 
         key= "Group",value = "TotAbundance") %>% 
  mutate(
    TimeGMT = hour(hms(TimeGMT)) + (minute(hms(TimeGMT)) + (second(hms(TimeGMT))/60))/60,
    TimeLocal = TimeGMT + round(Longitude/15), # Add 1hr for every 15 degrees of longitude
    TimeLocal = case_when(
      TimeLocal >= 0 & TimeLocal <= 24 ~ TimeLocal,
      TimeLocal > 24 ~ TimeLocal - 24,
      TimeLocal < 0 ~ TimeLocal + 24),
    Project = as.character(Project), 
    Type = case_when(
      Project == "cpr" | Project == "SO" ~ "CPR",
      Project == "nrs" ~ "Net"),
    Project = case_when(
      Project == "cpr" ~ "IMOS-CPR",
      Project == "nrs" ~ "IMOS-NRS",
      Project == "SO" ~ "SO-CPR"),
    Type = as.factor(Type), 
    Project = as.factor(Project),
    Group = as.factor(Group),
    DBase = as.factor("IMOS"),  # Assign level to databbase origin
    Transect = str_extract(ShipCruise,"[A-Za-z-]+")) %>%
  add_column(Volume=NA, Depth=NA, Upper_Z=NA, Lower_Z=NA, 
             Tow=NA, Mesh=NA, Ship=NA ) %>% # Add missing variables
  mutate(UniqueSampleID = str_c(ShipCruise, 1:n() ,sep="_"),
         Gear = case_when(Project == "IMOS-CPR" ~ 191,
                          Project == "SO-CPR" ~ 191,
                          Project == "IMOS-NRS" ~ 1000),
         Tow = case_when(Project == "IMOS-CPR" ~ "S",
                         Project == "SO-CPR" ~ "S",
                         Project == "IMOS-NRS" ~ "V")) %>% 
  dplyr::select(Retain)


SO <- imos %>% 
  filter(Project == "SO-CPR")

Other <- imos %>% 
  filter(Project != "SO-CPR") %>% 
  mutate(Transect = as.character(Transect))

# Get Southern Ocean Metadata
meta <- read.csv("Data/SO_CPR_RouteInfo_AAD.csv")
meta <- meta %>% 
  dplyr::select(c("TowNumber", "ShipCode"))

SO <- SO %>% mutate(
  ShipCruise = str_remove(ShipCruise,"_SO-CPR"),
  ShipCruise = as.numeric(ShipCruise))

SO2 <- left_join(SO,meta,by = c("ShipCruise" = "TowNumber"))
SO2 <- SO2 %>% 
  distinct() %>% 
  mutate(
    Transect = ShipCode,
    Transect = as.character(Transect),
    Transect = replace_na(Transect, "SOCPR_Temp"),
    ) %>% 
  dplyr::select(-ShipCode)

SO2 <- SO2 %>% mutate(
  ShipCruise = as.character(ShipCruise),
  ShipCruise = str_c(ShipCruise,"_SO-CPR"),
  ShipCruise = as.factor(ShipCruise))
                          

imos <- rbind(Other, SO2)
imos$Transect <- as.factor(imos$Transect)

# save data
saveRDS(imos,'IMOS_Final.rds')

rm(imos)
