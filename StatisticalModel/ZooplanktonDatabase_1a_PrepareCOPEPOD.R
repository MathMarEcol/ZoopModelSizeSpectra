rm(list = ls())

library(tidyverse)
library(lubridate)

source("fZerosByUniqueCruise.R")

Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")

######################### 1. READ RAW COPEPOD DATA ####################################
# The codes come from 
# https://www.st.nmfs.noaa.gov/copepod/2014/documents/fst37_copepod-2014.pdf
# The Plankton Grouping Code (PGC) is a seven digit number which identifies 
# the plankton taxa’s membership in up to four groups 
dat <- read_csv("Data/COPEPOD/copepod__4000000-compilation.csv", 
                na = c("null","----","-----","-99.000","n/a","-"))

dat <- dat %>% 
  mutate_all(list(~ifelse(.==-99, NA, .))) %>% 
  mutate_at(vars(TIMEgmt),list(~ifelse(.==99, NA, .))) %>% 
  mutate_at(vars(TIMEgmt),list(~ifelse(.==99.99, NA, .))) %>% 
  mutate_all(list(~ifelse(.==-999.9, NA, .))) %>%
  mutate(Group = NA,
         Depth = NA,
         NMFS_PGC = abs(NMFS_PGC)) %>% 
  rename(Longitude = LONGITDE, Latitude = LATITUDE, 
         ShipCruise = "SHP-CRUISE", Year = YEAR, Month = MON, 
         Day = DAY, TimeGMT = TIMEgmt, TimeLoc = TIMEloc, Upper_Z = UPPER_Z, 
         Lower_Z = LOWER_Z, Tow = T, Gear = GEAR, Mesh = MESH, 
         Volume = "Water Strained", OrigValue = "Original-VALUE", OrigUnits = "Orig-UNITS",
         Abundance = "VALUE-per-volu", Units = UNITS, SpeciesName = "SCIENTIFIC NAME -[ modifiers ]-", 
         Record.ID = "RECORD-ID", Project = "DATASET-ID", Ship = SHIP, Inst = INST) %>% 
  mutate(Transect = as.factor(ShipCruise)) %>% 
  filter(V == "c") %>% 
  filter(Gear != 222) %>% 
  droplevels()

# 
# dat <- readRDS("LatestDatabaseOuput_Final_Enviro.rds") # n=980,90
# # datx <- dat %>% filter(TotAbundance > 100000)
# datx <- dat %>% filter(TimeLocal == TimeGMT)
# datx <- datx %>% filter(abs(Longitude)>15)
# unique(datx$Project)
# datx <- dat %>% filter(Project == "North Pacific and Bering Sea Oceanography 1958")
# 
# # Assume that time is recorded in local time
# plot(datx$TimeLocal,datx$TotAbundance)
# 
# # Assume that time is actually in GMT
# datx$TimeLocal <- datx$TimeLocal + round(datx$Longitude/15)
# datx <- datx %>%   mutate(TimeLocal = case_when(TimeLocal >= 0 & TimeLocal <= 24 ~ TimeLocal,
#                                                 TimeLocal > 24 ~ TimeLocal - 24,
#                                                 TimeLocal < 0 ~ TimeLocal + 24))
# plot(datx$TimeLocal,datx$TotAbundance)


## I think they are GMT so I need to change to local time. Set GMT to NA so it's caught below
dat <- dat %>% 
  mutate(TimeGMT = replace(TimeGMT, Project == "North Pacific and Bering Sea Oceanography 1958", NA)) %>% 
  mutate(TimeGMT = replace(TimeGMT, TimeGMT <0 | TimeGMT >24, NA)) %>%
  mutate(TimeLocal = case_when(TimeLoc >= 0 ~ TimeLoc,
                               TimeGMT >= 0 & TimeGMT <= 24 & is.na(TimeLoc) ~ TimeGMT + round(Longitude/15),
                               is.na(TimeLoc) & is.na(TimeGMT) ~ runif(is.na(TimeLoc) & is.na(TimeGMT))*24)) %>% 
  mutate(TimeLocal = case_when(TimeLocal >= 0 & TimeLocal <= 24 ~ TimeLocal,
                               TimeLocal > 24 ~ TimeLocal - 24,
                               TimeLocal < 0 ~ TimeLocal + 24)) %>%
  mutate_at(vars(Mesh),list(~ifelse(Project == "HAKUHO MARU collection",100, .))) %>% #https://watermark.silverchair.com/18-5-673.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAkQwggJABgkqhkiG9w0BBwagggIxMIICLQIBADCCAiYGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMnj91yrSRXnnzZx2sAgEQgIIB98nTygHOEl2a-qVSn2V4pQDE5ITo62SqD4ChmeDO9bzAJRwI-RdlpTY17VVbLXYYSacw1sE1TWqb7P0Sab_7FmUwSXLa950aubx7dbAIla2c9v6k9qNMXXt0Ce3d3QgJeH_MiYdTAleWkdQikv_o7YTM3aagEYP5pEotdB9UrQ4y2fwlHRXAxoGesZq-eZ_iEE_AY2q3CP1EzdeeNaQfozfh9ZcLqTZNYxGXpQ42Z_81Q_3VThvwRuywOalICSlapxbOW6-R374EKBbZuXORjcZUZwEKVuxFjFDXBD35wBgBTT5IcbnTptb9Sbw8rTcMAwKy3HZYHTig5PpbQ6LWDRtRa1aePS7chFbORWwfZiGbuLuA_Sl0UASN3wPrWlusKFq0-Wcr_mgdqlCA1bQNt5xpIWwUB9qcnraxEvMDgE0G-6EZGrmYBiOPCXSDEUZSfmNkJEaLrJUJTEA8zmEJDftJaYXP7Kj-HREAaJyjF8qaqqp6G5-SERz74ZqeBHku3Scvn1ur0c7fqRk0PvIkVIq_XGY79-ar-QIsJtLESqQvjXVnNWS_DzdQEu37641QUPA3YYqs_Nywk088bXCQ8IwP0rOT_Oys5k9y7nndOGTh3vBwqBqKC1QivN_LTb814EYbQLgiqQH1L571SHl-03x2DfcNBcg4
  mutate_at(vars(Group),list(~ifelse(NMFS_PGC >= 4212000 & NMFS_PGC < 4213000, "Copepod", .))) %>% 
  mutate_at(vars(Group),list(~ifelse(NMFS_PGC == 4320000, "Chaetognaths", .))) %>% 
  mutate_at(vars(Group),list(~ifelse(NMFS_PGC == 4218000, "Euphausiids", .))) %>% 
  mutate_at(vars(Group),list(~ifelse(NMFS_PGC >= 4032000 & NMFS_PGC < 4038000, "Jellyfish", .))) %>% 
  mutate_at(vars(Group),list(~ifelse(NMFS_PGC == 4357500, "Larvaceans", .))) %>% 
  mutate_at(vars(Group),list(~ifelse(NMFS_PGC == 4355010, "Salps", .))) %>% 
  subset(!is.na(Group)) %>% 
  filter(!is.na(Abundance)) %>% # mutate(Units = as.factor(Units))
  filter(!is.na(Units) & !is.na(TimeLocal)) %>% 
  filter(Project != "SAHFOS-CPR Atlantic Ocean" & Project != "EcoMon-SOOP (Mid-Atlantic Bight)" & Project != "EcoMon-SOOP (Gulf of Maine)") %>% 
  filter(Gear != 127) %>% # Remove bottle samples
  dplyr::select(-c("F1","F2","F3","F4","F1_1","F2_1","F3_1","F4_1","X42")) %>%  # Neaten things up
  mutate(UniqueSampleID = paste0(ShipCruise, Record.ID), # Create Unique ID for every net deployment
         UniqueSampleID = str_sub(UniqueSampleID, end = -6), # Remove last 5 characters which are the 'taxaID')
         UniqueSampleID = paste0(UniqueSampleID, sprintf('%04.0f', Upper_Z)), # Now add in depth as well
         UniqueSampleID = as.factor(UniqueSampleID)) %>% # Unique ID for each sample
  arrange(Group, UniqueSampleID) %>%  # Arrange so by Group and UniqueSampleID
  filter(str_detect(OrigUnits, "#/m3") | str_detect(OrigUnits, "#/haul")) %>% 
  mutate(DBase = as.factor("COPEPOD"), # Assign level to databbase origin
         Type = "Net", # All COPEPOD samples are from net/trawls. We remove CPR
         ShipCruise = as.factor(ShipCruise),
         Project = as.factor(Project),
         Group = as.factor(Group),
         Volume = as.numeric(str_remove(Volume," m3"))) # Remove m3 from the column

datx <- dat %>% 
  group_by(UniqueSampleID) %>% 
  summarise(sd = sd(Volume,na.rm = T), # Find all samples with more than 1 volume
            Project = Project[1]) %>% 
  filter(sd>0) %>% 
  ungroup() %>% 
  distinct(UniqueSampleID)

dat <- dat %>% 
  filter(!UniqueSampleID %in% datx$UniqueSampleID) # Remove all the sampleIDs which have dodgy volume recordings
rm(datx)

######## 2. Create Omnivorous and Carnivorous copepod groups ###############
# Import in FeedingType and clean it
FeedingType <- read_csv("Data/Copepod_Feeding_Type_v3_wAPHIA.csv")
FeedingType <- FeedingType %>% 
  dplyr::select(-c(TAXON_NAME, aphiaID)) %>%
  distinct()

dat <- dat %>% 
  filter(Group == "Copepod") %>% 
  mutate(Genus = str_extract(SpeciesName,"^\\S*"),
         Genus = str_replace(Genus, "Microcalanids", "Microcalanus"), # MisSpelling
         Genus = str_replace(Genus, "Scolecithrichopsis", "Scolecitrichopsis"), # MisSpelling
         Genus = str_replace(Genus, "Undinopsis","Bradyidius"), # New name
         Genus = str_replace(Genus, "Idyaea","Tisbe"), # New name
         Genus = str_replace(Genus, "Diaxis","Diaixis"), # MisSpelling            
         Genus = str_replace(Genus, "P-cal", "Paracalanus")) %>% 
  left_join(FeedingType, by = c("Genus" = "GENUS_NAME")) %>% 
  mutate(Diet = case_when(FEED == "CC" ~ "Carnivore", # Diet: Carnivores (CC) vs 
                          FEED == "CBF" ~ "Omnivore", # Omnivores (CBF + CD + CH + CO + CP + SF)
                          FEED == "CD" ~ "Omnivore",
                          FEED == "CH" ~ "Omnivore",
                          FEED == "CO" ~ "Omnivore",
                          FEED == "CP" ~ "Omnivore",
                          FEED == "SF" ~ "Omnivore"),
         Diet = as.factor(Diet),
         Group = Diet,
         Group = str_replace(Group, "Omnivore","OmniCopepods"),
         Group = str_replace(Group, "Carnivore","CarnCopepods")) %>% 
  filter(!is.na(Group)) %>% 
  mutate(Diet2 = case_when(FEED == "CC" ~ "Carnivore", # Diet2: vs Carnivore (CC + CBF + CD + CH + CO + CP + SF)
                           FEED == "CBF" ~ "Carnivore",
                           FEED == "CD" ~ "Carnivore",
                           FEED == "CH" ~ "Carnivore",
                           FEED == "CO" ~ "Carnivore",
                           FEED == "CP" ~ "Carnivore",
                           FEED == "SF" ~ "Carnivore",
                           FEED == "CH" ~ "Herbivore"), # Herbivore (CH) 
         Diet2 = as.factor(Diet2)) %>% 
  dplyr::select(-c("Diet", "Diet2","FEED","Genus")) %>% 
  rbind(dat) %>%  # Merge grouped copepods back in
  filter(Group != "Copepod") %>% # Remove copepods from dat, 
  droplevels() %>% 
  mutate(Group = as.factor(Group),
         UniqueSampleID = as.factor(UniqueSampleID))%>%
  group_by(Group, UniqueSampleID)


######## 3 Merge copepods and create totals ###############
dat_sum <- dat %>% 
  dplyr::summarise(TotAbundance = sum(as.numeric(Abundance)))# Calculate total abundance for each UniqueSampleID

# All rows in a UniqueSampleID are the same except for ScientificName, so keep the distinct ones 
dat <- dat %>% distinct(Group, UniqueSampleID, .keep_all = TRUE) %>% 
  right_join(dat_sum, by = c("Group", "UniqueSampleID")) %>%
  ungroup() %>% 
  mutate(Group = as.factor(Group),
         UniqueSampleID = as.factor(UniqueSampleID)) %>% 
  dplyr::select(Retain)

######## 4 Add Zeros ###############
uni <- unique(dat$ShipCruise) # Get unique cruises

dat0 <- fZerosByUniqueCruise(dat,uni)


# save data
saveRDS(dat0,'COPEPOD_Final.rds')

rm(dat, dat_sum, FeedingType, Retain)
