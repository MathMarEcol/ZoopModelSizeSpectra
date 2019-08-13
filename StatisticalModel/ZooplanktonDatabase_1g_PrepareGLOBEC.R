# "Zooplankton were collected using a 0.5 meter diameter, 202 micron mesh net towed 
# verticaly from within 5 meters of the bottom (to a maximum depth of 100 m) to the 
# surface at a rate of 30 meters per minute. A TSK flowmeter was used to monitor the 
# amount of water filtered. The samples were preserved in a 5% buffered formalin/seawater 
# solution." (W.T. Peterson, et al., 2002 p. 392)
# "In the laboratory, the zooplankton samples were diluted and subsampled with a 1.1 ml 
# Stempel pipette. Two to four such subsamples, about 2% of the total sample, were counted 
# at 25-50x magnification. Copepods and Euphausiids were identified to species and 
# developmental stage; other zooplankton were assigned to broad taxonomic groups (e.g. 
# polychaetes, medusae, larvaceans, chaetognaths.) All euphausiids, pteropods, salps, and 
# chaetognaths were measured. In each sample, the population density of each taxonomic 
# group (number of individuals per cubic meter) was calculated. Copepod densities were 
# converted to biomass estimates using dry weight/developmental-stage values found in the 
# literature; biomasses of euphausiids, pteropods, salps, and chaetognaths were calculated 
# from densities using length-weight regressions found in the literature." (W.T.Peterson, 
# et al., 2002, p.392)

source("func/fZerosByUniqueCruise.R")

library(ncdf4)
library(tidyverse)
library(lubridate)

Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")

nc <- nc_open("~/Dropbox/EcosystemEfficiency/StatisticalModel/ZoopAbundanceModels/Data/GLOBEC/Globec_vpt.nc")
Latitude = as.numeric(ncvar_get(nc, "lat"))
dat <- as.data.frame(Latitude)

dat <- dat %>% 
  add_column(Longitude = as.numeric(ncvar_get(nc, "long")),
             ShipCruise = as.factor(ncvar_get(nc, "cruise_id")),
             Month = as.numeric(ncvar_get(nc, "month_local")),
             Day = as.numeric(ncvar_get(nc, "day_local")),
             Year = as.numeric(ncvar_get(nc, "year")),
             Gear = 330, # VPT Vertical Plankton Tow.
             Mesh = as.numeric(ncvar_get(nc, "gear_mesh")) * 1e3,
             Abundance = as.numeric(ncvar_get(nc, "abund_m3")),
             Species = ncvar_get(nc, "genus_species"),
             Project = as.factor(ncvar_get(nc, "program")),
             Volume = as.numeric(ncvar_get(nc, "vol_filt")),
             Upper_Z = as.numeric(ncvar_get(nc, "min_sample_depth")),
             Lower_Z = as.numeric(ncvar_get(nc, "max_sample_depth"))) %>% 
  mutate(TimeRaw = hms(ncvar_get(nc, "time_local_time")),
         TimeLocal = hour(TimeRaw) + (minute(TimeRaw) + second(TimeRaw)/60)/60,
         Depth = (Lower_Z-Upper_Z)/2,
         Type = "Net",
         Ship = ShipCruise,
         Transect = ShipCruise,
         Tow = "O")


Feed <- read_csv("Data/Copepod_Feeding_Type_v3_wAPHIA.csv")
jellies = str_to_upper(c("Muggiaea","Leuckartiara","Diphyidae","Anthomedusae","Phialidium","Tiaropsis","Medusa",
                         "Obelia","Aequorea","Leptomedusae","Mitrocoma","Eutonina","Praya","Halitholus", "Hydromedusa"))

dat <- dat %>% 
  add_column(DBase = "GLOBEC", TimeGMT = NA) %>% 
  mutate(Taxa = str_to_title(str_extract(dat$Species,"^[^_]+"))) %>% 
  left_join(Feed, by=c("Taxa" = "GENUS_NAME")) %>% 
  mutate(Group = case_when(FEED == "CC" ~ "CarnCopepods", # Diet: Carnivores (CC) vs 
                        FEED == "CBF" ~ "OmniCopepods", # Omnivores (CBF + CD + CH + CO + CP + SF)
                        FEED == "CD" ~ "OmniCopepods",
                        FEED == "CH" ~ "OmniCopepods",
                        FEED == "CO" ~ "OmniCopepods",
                        FEED == "CP" ~ "OmniCopepods",
                        FEED == "SF" ~ "OmniCopepods",
                        str_detect(Species,"SALP") ~ "Salps",
                        str_detect(Species,"THETYS_VAGINA") ~ "Salps",
                        str_detect(Species,"LARVACEA") ~ "Larvaceans",
                        str_detect(Species,"EUPHAUS") ~ "Euphausiids",
                        str_detect(Species,"THYSANOESSA") ~ "Euphausiids",
                        str_detect(Species,"STYLOCHEIRON") ~ "Euphausiids",
                        str_detect(Species,"NYCTIPHANES") ~ "Euphausiids",
                        str_detect(Species,"TESSARABRACHION") ~ "Euphausiids",
                        str_detect(Species,"NEMATOBRACHION") ~ "Euphausiids",
                        str_detect(Species,"NEMATOSCELIS") ~ "Euphausiids",
                        str_detect(Species,"CHAETOGNATHA") ~ "Chaetognaths",
                        str_detect(Species,"JELLYFISHES") ~ "Jellyfish",
                        str_detect(Species,paste(jellies,collapse = '|')) ~ "Jellyfish"))

dat <- dat %>% 
  select(-c("aphiaID", "FEED", "TAXON_NAME", "Species", "TimeRaw", "Taxa")) %>% 
  distinct(.keep_all = FALSE) %>% 
  filter(!is.na(Group))
  
dat <- dat %>% 
  mutate(UniqueSampleID = paste0(ShipCruise,Year,str_pad(Month, 2, pad = "0"),Day,sprintf('%02.4f', TimeLocal)),
         UniqueSampleID = str_replace(UniqueSampleID,"\\.",""),
         UniqueSampleID = str_sub(UniqueSampleID, end = -2),
         UniqueSampleID = as.factor(UniqueSampleID),
         Group = as.factor(Group)) %>% 
  group_by(UniqueSampleID, Group)

  
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


## Add zeros to the data by cruise
uni <- unique(dat$ShipCruise) # Get unique cruises

dat0 <- fZerosByUniqueCruise(dat,uni)

dat <- dat0 %>% 
  select(Retain)

# save data
saveRDS(dat,'DatabaseOutput/GLOBEC_Final.rds')


# THE CODE BELOW IS TO FIND THE CLASSIFICATION INFO FOR ALL THE UNMATCHED TAXA
# I HAVE ALREADY CHECKED AND THER EIS NOTHING RELAVENT FOR US IN THIS LIST
# SO NO NEED TO USE THIS CODE AGAIN.
# 
# # Get the rows which don't match
# temp <- dat0$Taxa[is.na(dat0$FEED)]
# temp <- unique(temp)
# 
# library(taxize)
# library(worrms)
# Class=list() 
# 
# for(i in 1:length(temp)) {
#   tryCatch(
#     {
#       aphiaID <- wm_records_taxamatch(temp[i], marine_only = FALSE, failonerror = FALSE)  
#       Class[i] <- aphiaID[[1]]$class
#       rm(aphiaID)
#       },
#   error=function(e) {
#     Class[i] <- NA
#   })
# }
# 
# 
# Class[sapply(Class, is.null)] <- NA
# Class2 <- unlist(Class)
# 
# check <- as.data.frame(Class2)
# check$Taxa <- temp
