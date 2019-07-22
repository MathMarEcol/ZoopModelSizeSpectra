library(tidyverse)
library(lubridate)

Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")

salps <- read_csv("Data/Henschke/SalpDatabase2008-2010.csv")

salps <- salps %>% 
  rename(TimeGMT = `Time (UTC)`, SampleID = `CTD station`,
         TotAbundance = `Numeric density (ind. m-3)`, Lower_Z = `Sampling Depth (m)`,
         Transect = Voyage) %>% 
  select(-c(Phylum, Class, Order, Family, Genus, Species, Taxon, `Presence/Absence`, 
            `Water Depth (m)`, `Net opening (cm)`, `Collection method`)) %>% 
  add_column(Tow = "V", Type = "Net", Ship = as.factor("SouthernSurveyor"), DBase = "Henschke", 
             Depth = 25, Mesh = 225, Group = as.factor("Salps"), Upper_Z = 0, Volume = 19.24226) %>% 
  mutate(TimeGMT = hour(hms(TimeGMT)) + (minute(hms(TimeGMT)) + (second(hms(TimeGMT))/60))/60,
         TimeLocal = TimeGMT + round(Longitude/15),
         TimeLocal = case_when(TimeLocal >= 0 & TimeLocal <= 24 ~ TimeLocal,
                               TimeLocal > 24 ~ TimeLocal - 24,
                               TimeLocal < 0 ~ TimeLocal + 24),
         ShipCruise = as.factor(str_c(Ship,as.character(Year),sep = "_")),
         Transect = as.factor(Transect),
         Project = Transect,
         Gear = 135) %>%  # 135 is the COPEPOD gear code for N70 
  mutate(Date = dmy(Date),
       Day = day(Date), 
       Month = month(Date)) %>% 
  select(-c(Date))

salps <- salps %>% 
  mutate(SampleID = as.character(SampleID),
         SampleID = str_pad(SampleID, 2, pad = "0"),
         SampleID = str_c(as.character(Year),SampleID,sep = "_")) %>% 
  group_by(SampleID) %>% 
  mutate(id = row_number(),
         id = as.character(id),
         UniqueSampleID = str_c(SampleID,id,sep = "_")) %>% 
  ungroup() %>% 
  select(-SampleID)

# All rows in a IDforSumming are the same except for ScientificName, so keep the distinct ones
salps <- salps %>% mutate(Group = as.factor(Group)) %>%
  dplyr::select(Retain) %>% 
  droplevels()

saveRDS(salps,'Henschke_Final.rds')
