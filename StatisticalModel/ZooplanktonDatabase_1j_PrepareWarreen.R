library(tidyverse)
library(lubridate)

Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")

warr <- read_csv("Data/Warreen/warreen4csiro.csv", na = "NaN")

warr <- warr %>% 
  rename(Day = day, Month = month, Year = year, Latitude = latitude, 
         Longitude = longitude, Station = `St no.`, TimeLocalx = `Time-of-day`,
         Tow = `Net type`, TowTime = `Tow time [min]`,
         Euphausiids = euphausiids, Salps = thaliacea, 
         Larvaceans = larvaceans, Chaetognaths = chaetognaths) %>% 
  select(-c(`Cloud-free solar rad [W/m2]`, SST, `CARS SST`, SSS, 
            `CARS SSS`, `20th Cent wind u [m/s]`, `20th Cent wind v [m/s]`,
            `Biovolume [cc]`, `Dominant organism`)) %>% 
  add_column(Type = "Net", Ship = as.factor("Warreen"), DBase = "Warreen", 
             Mesh = 222, Upper_Z = 0, TimeGMT = NA, Project = as.factor("Warreen"),
             Gear = 135, Transect = as.factor("Warreen"), ShipCruise = as.factor("Warreen")) %>%  # 135 is the COPEPOD gear code for N70 
  mutate(Tow = case_when(Tow == "N70_V50_0" ~ "V",
                         Tow == "N70_H_0" ~ "H"),
         Lower_Z = case_when(Tow == "V" ~ 50,
                             Tow == "H" ~ 0),
         Depth = case_when(Tow == "V" ~ 25,
                             Tow == "H" ~ 0),
         Volume = case_when(Tow == "V" ~  pi * 0.35 *0.35 * 50, # Area of net and depth of tow (50 m)
                            Tow == "H" ~ TowTime * 60 * 0.89 * pi * 0.35 *0.35), # Based on net area, approx speed and approx time (Accurate to about 15 minutes)
         UniqueSampleID = str_c(as.character(Year),Tow,str_pad(Station, 3, pad = "0"),sep = "_"),
         ShipCruise = as.factor(str_c(Ship,as.character(Year),sep = "_"))) %>% 
  filter(Tow == "V" | Tow == "H")

warr <- warr %>% 
  mutate(TimeLocal = as.numeric(TimeLocalx),
         TimeLocal = replace(TimeLocal, TimeLocal>2400, NA))

warr <- warr %>% 
  mutate(TimeLocal = replace(TimeLocal, is.na(TimeLocal), runif(sum(is.na(TimeLocal)), min = 0, max = 24)*100),
         H = floor(TimeLocal/100),
         M = TimeLocal-H*100,
         TimeLocal = H + M/60) %>% 
  select(-TimeLocalx)




warr <- warr %>% 
  gather(c("Chaetognaths", "Euphausiids", "Larvaceans", "Salps"), 
         key= "Group",value = "TotAbundance") %>% 
  mutate(TotAbundance = TotAbundance/(Volume/10))
         
warr <- warr %>% mutate(Group = as.factor(Group)) %>%
  dplyr::select(Retain) %>% 
  droplevels()

saveRDS(warr,'Warreen_Final.rds')
