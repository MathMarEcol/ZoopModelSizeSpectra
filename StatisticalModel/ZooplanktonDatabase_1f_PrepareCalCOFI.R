rm(list = ls())
source("func/fZerosByUniqueCruise.R")

library(tidyverse)
library(lubridate)

Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")

dat1d <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_Carnivores_Day.csv",skip=6, na = c("N/A", "NA"))
dat1d <- dat1d %>% mutate(
  Group = "OmniCopepods",
  DN = "Day") %>% 
  select(-c("Tows","Source"))

dat1n <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_Carnivores_Night.csv",skip=6, na = c("N/A", "NA"))
dat1n <- dat1n %>% mutate(
  Group = "OmniCopepods",
  DN = "Night") %>% 
  select(-c("Tows","Source"))

dat2d <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_Omnivores_Day.csv",skip=6, na = c("N/A", "NA"))
dat2d <- dat2d %>% mutate(
  Group = "CarnCopepods",
  DN = "Day") %>% 
  select(-c("Tows","Source"))

dat2n <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_Omnivores_Night.csv",skip=6, na = c("N/A", "NA"))
dat2n <- dat2n %>% mutate(
  Group = "CarnCopepods",
  DN = "Night") %>% 
  select(-c("Tows","Source"))

dat3d <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_APPENDICULARIA_Day.csv",skip=6, na = c("N/A", "NA"))
dat3d <- dat3d %>% mutate(
  Group = "Larvaceans",
  DN = "Day") %>% 
  select(-c("Tows","Source"))

dat3n <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_APPENDICULARIA_Night.csv",skip=6, na = c("N/A", "NA"))
dat3n <- dat3n %>% mutate(
  Group = "Larvaceans",
  DN = "Night") %>% 
  select(-c("Tows","Source"))

dat4d <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_Euphausiids_Day.csv",skip=6, na = c("N/A", "NA"))
dat4d <- dat4d %>% 
  mutate(Group = "Euphausiids",
         DN = "Day")

dat4n <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_Euphausiids_Night.csv",skip=6, na = c("N/A", "NA"))
dat4n <- dat4n %>% 
  mutate(Group = "Euphausiids",
         DN = "Night")

dat5d <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_SALPIDA_Day.csv",skip=6, na = c("N/A", "NA"))
dat5d <- dat5d %>% mutate(
  Group = "Salps",
  DN = "Day") %>% 
  select(-c("Tows","Source"))

dat5n <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_SALPIDA_Night.csv",skip=6, na = c("N/A", "NA"))
dat5n <- dat5n %>% mutate(
  Group = "Salps",
  DN = "Night") %>% 
  select(-c("Tows","Source"))

dat6d <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_CHAETOGNATHA_Day.csv",skip=6, na = c("N/A", "NA"))
dat6d <- dat6d %>% mutate(
  Group = "Chaetognaths",
  DN = "Day") %>% 
  select(-c("Tows","Source"))

dat6n <- read_csv("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZooplanktonModels/Data/CalCOFI/CalCOFI_CHAETOGNATHA_Night.csv",skip=6, na = c("N/A", "NA"))
dat6n <- dat6n %>% mutate(
  Group = "Chaetognaths",
  DN = "Night") %>% 
  select(-c("Tows","Source"))


dat <- rbind(dat1d, dat1n, dat2d, dat2n, dat3d, dat3n, dat4d, dat4n, dat5d, dat5n, dat6d, dat6n)
# rm(dat1d, dat1n, dat2d, dat2n, dat3d, dat3n, dat4d, dat4n, dat5d, dat5n, dat6d, dat6n)

dat <- dat %>% 
  add_column(Volume=NA, Depth=NA, Upper_Z=0, Lower_Z=NA, 
             Tow="O", Mesh=505, TimeLocal=NA, TimeGMT = NA,
             DBase="CalCOFI", Type="Net", Project = "CalCOFI",
             Gear=NA) %>% 
  rename("Transect" = Line, "ShipCruise" = Cruise)

dat <- dat %>%
  filter(!is.na(Date)) %>% 
  mutate(Lower_Z = case_when(year(Date) <= 1950 ~ 70,
                             year(Date) > 1950 & year(Date) <= 1968 ~ 140,
                             year(Date) >= 1969 ~ 210),
         Gear = case_when(year(Date) <= 1977 ~ 274,
                          year(Date) >= 1978 ~ 275),
         Abundance = Abundance/1e3,
         Depth = (Lower_Z-Upper_Z)/2,
         Day = day(Date),
         Month = month(Date),
         Year = year(Date),
         Transect = as.factor(Transect),
         UniqueSampleID = paste0(ShipCruise,Station))

# A 1-m diameter ring net was used from 1949 - late 1977, which was then 
# replaced with a 0.71-m diameter bongo net in December 1977 (Ohman and Smith 1995). 
# Gear Types from COPEPOD
# 274	  CalCOFI standard 1-meter net (1951-1978)
# 275	  CalCOFI standard bongo net (1978-present)

      
# The maximum zooplankton sampling depth in the CalCOFI time series increased from 
# 0-70 m (1949-1950) to 
# 0-140 m in 1951. It was maintained at 0-140 m from 1951-1968, and then increased to 
# 0-210 m in 1969, where it has since remained (see Ohman and Smith 1995). 

# CalCOFI data is split into D/N only. Here we randomly assign a time.
dat <- dat %>%
  mutate(TimeLocal = case_when(DN=="Day" ~ (runif((DN=="Day"))*12)+6,
                               DN=="Night" ~ (runif((DN=="Night"))*-12)+30),
         TimeLocal = case_when(TimeLocal > 24 ~ TimeLocal-24,
                               TimeLocal < 24 ~ TimeLocal))


dat <- dat %>% group_by(UniqueSampleID, Group)


datsum <- dat %>% 
  dplyr::summarise(TotAbundance = sum(Abundance))

dat2 <- dat %>% 
  distinct(Group, UniqueSampleID, .keep_all = TRUE) %>% 
  right_join(datsum, by = c("Group", "UniqueSampleID")) %>%
  ungroup() %>% 
  mutate(Group = as.factor(Group),
         UniqueSampleID = as.factor(UniqueSampleID))

dat <- dat2
rm(datsum, dat2)


## Add zeros
uni <- unique(dat$ShipCruise)
dat0 <- fZerosByUniqueCruise(dat,uni)


dat <- dat0 %>%
  select(Retain)



# save data
saveRDS(dat,'DatabaseOutput/CalCOFI_Final.rds')





