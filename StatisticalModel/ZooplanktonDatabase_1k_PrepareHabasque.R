library(tidyverse)
library(lubridate)

# To help multifrequency acoustic data validation, zooplankton samples were 
# acquired near each PIRATA buoy during FR26 survey. Oblique tows were
# conducted with bongo nets (300 µm mesh), towed between on average 200 m 
# depth and the sea surface. A flowmeter (KC Denmark 23.090) mounted in
# the mouth of the net was used to measure seawater filtered during the tow. 
# Samples were fixed in 4% formaldehyde buffered with borax for laboratory 
# analysis and stored at -20°C. Mean filtered water volume was 175m3

Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")

habas <- read_delim("Data/Habasque/Habasque61651.csv", ";")

# I am guess that "inf" (inferior) is a FR abbrev for less than and "sup" (superior) is the abbrev for greater than
habas <- habas %>% 
  select(-starts_with("Fraction_")) %>% 
  mutate_at(vars(starts_with("Counts")), list(~ifelse(is.na(.), 0, .))) %>% 
  mutate(Abundance = rowSums(select(., starts_with("Counts")))) %>% 
  select(-starts_with("Counts"))

habas <- habas %>% 
  mutate(Date = dmy(Date),
         Day = day(Date),
         Month = month(Date),
         Year = year(Date))

habas <- habas %>% 
  mutate(Time = str_replace(Time,"h$","h00"), # If h is last, it means not minsutes so lets add them
         Time2 = hm(Time),
         TimeLocal = hour(Time2) + (minute(Time2)/60)) %>% 
  select(-c(Time,Time2, Group))

habas <- habas %>% 
  mutate(Genus = str_extract(Taxa,"^\\S*"),
         Genus = str_replace(Genus, "Microcalanids", "Microcalanus"), # MisSpelling
         Genus = str_replace(Genus, "Saphirinia", "Sapphirina"), # MisSpelling
         Genus = str_replace(Genus, "Clausocalanus/Paracalanus", "Clausocalanus"), # MisSpelling
         Genus = str_replace(Genus, "Undinulla", "Undinula"), # MisSpelling
         Genus = str_replace(Genus, "Sub/Eucalanus", "Eucalanus"), # MisSpelling
         Genus = str_replace(Genus, "Centrogaptilus", "Centraugaptilus"), # MisSpelling
         Genus = str_replace(Genus, "Chiridus", "Chiridius")) # MisSpelling

# Get the copepods seperately so we can check feeding pref
copepod <- habas %>% 
  filter(Order == "Copepod")

# Store everything else
other <- habas %>%
  filter(Order != "Copepod") %>% 
  mutate(Group = case_when(str_detect(Order,"Euphausiid") ~ "Euphausiids",
                           str_detect(Order,"Chaetognaths") ~ "Chaetognaths",
                           str_detect(Order,"Appendicularians") ~ "Larvaceans",
                           str_detect(Order,"Jellyfish") ~ "Jellyfish",
                           str_detect(Order,"Salps") ~ "Salps")) %>% 
  filter(!is.na(Group)) %>% 
  select(-c(Taxa, Genus, Order))

rm(habas)

FeedingType <- read_csv("Data/Copepod_Feeding_Type_v3_wAPHIA.csv")
FeedingType <- FeedingType %>% 
  dplyr::select(-c(TAXON_NAME, aphiaID)) %>%
  distinct()

copepod <- copepod %>% 
  left_join(FeedingType, by = c("Genus" = "GENUS_NAME")) %>% 
  mutate(Group = case_when(FEED == "CC" ~ "CarnCopepods", # Diet: Carnivores (CC) vs 
                           FEED == "CBF" ~ "OmniCopepods", # Omnivores (CBF + CD + CH + CO + CP + SF)
                           FEED == "CD" ~ "OmniCopepods",
                           FEED == "CH" ~ "OmniCopepods",
                           FEED == "CO" ~ "OmniCopepods",
                           FEED == "CP" ~ "OmniCopepods",
                           FEED == "SF" ~ "OmniCopepods"),
         Group = as.factor(Group)) %>% 
  filter(!is.na(Group)) %>% 
  select(-c(FEED, Taxa, Genus, Order))

habas <- rbind(copepod,other)

habas <- habas %>% 
  rename(ShipCruise = Survey, Lower_Z = Depth_max) %>% 
  add_column(Type = "Net", Ship = as.factor("Thalassa"), DBase = "Habasque", 
             Mesh = 00, Upper_Z = 0, TimeGMT = NA, Gear = 118, Tow = "O", Volume = 175) %>%  # 118 = Bongo in COPEPOD
  mutate(Transect = as.factor(ShipCruise),
         Project = ShipCruise,
         ShipCruise = as.factor(str_replace(ShipCruise,"PIRATA FR26", "PIRATA_FR26")),
         UniqueSampleID = as.factor(str_c(ShipCruise,Station,sep = "_")),
         Depth = (Lower_Z-Upper_Z)/2,
         Abundance = Abundance / Volume)

habas_sum <- habas %>% 
  group_by(UniqueSampleID, Group) %>% 
  dplyr::summarise(TotAbundance = sum(as.numeric(Abundance))) # Calculate total abundance for each UniqueSampleID

habas <- habas %>% 
  select(-Abundance)

habas <- left_join(habas, habas_sum) %>% 
  distinct()

rm(habas_sum)

habas <- habas %>% 
  mutate(Group = as.factor(Group)) %>%
  dplyr::select(Retain) %>% 
  droplevels()

saveRDS(habas,'DatabaseOutput/Habasque_Final.rds')


