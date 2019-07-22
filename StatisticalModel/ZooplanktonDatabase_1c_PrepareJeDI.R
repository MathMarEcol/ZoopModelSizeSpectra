library(yaImpute)
library(tidyverse)
library(lubridate)


Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")


jedi <- read_csv("Data/JeDI/JEDI_v3.csv",na="-9999.0",col_types = "ccccccdddtddcccccccffddddddddddddfffd")


jedi <- jedi %>% 
  filter(!is.na(date)) %>%
  # filter(!is.na(time_local) %>%
  mutate(date = mdy(date)) %>%
  mutate_each(list(~replace(., . == -9999, NA)), -project_title) %>% # project_title is here only because I couldn't get it to run for all columns
  rename(Year = year, Month = month, Day = day, TimeLocal = time_local,
         Latitude = lat, Longitude = lon,
         Depth = depth, Upper_Z = depth_upper,
         Lower_Z = depth_lower, Gear = collection_method, Mesh = net_mesh,
         Abundance = density,  Project = owner_dataset,
         Inst = project_title, data_type=data_type) %>%
  mutate(DBase = as.factor("JeDI"), # Assign level to databbase origin
         Mesh = Mesh*1e3, # Mesh is in mm. Convert to um
         ShipCruise = Inst, # Use project as the Ship-Cruise random variable
         UniqueSampleID = Inst, # Use project as the UniqueSampleID random variable
         Transect = as.factor(ShipCruise)) %>% 
  add_column(Volume=NA, OrigValue=NA, OrigUnits="m3", Units="m3", # Add missing variables
             Tow=NA, TimeGMT=NA, Type="Net", Record.ID=NA, Ship=NA)

# Filter the data to keep what I want
jedi <- jedi %>% 
  filter(Project != "Lucas_C" & #Lucas data is within a lake.
           ShipCruise != "L4" & # L4 processed seperately and  Remove.
           (Gear == "plankton_net" | Gear == "trawl") & # Only net and trawl surveys)
           data_type == "quantitative" & # Only quantitative surveys
           Inst != 'COPEPOD' & # No surveys from the COPEPOD database
           # !is.na(Mesh) & # Only surveys with mesh size recorded
           !is.na(Abundance) & # Filter for missing abundance
           # !is.na(TimeLocal) & # Filter for missing Time # There actually aren't any: runif(is.na(TimeLoc) & is.na(TimeGMT))*24))
           (taxon == "salp" | rank_phylum == "Cnidaria" | taxon == "medusa") &
           taxon != "siphonophore" &
           Inst != "VIMS_Survey" & # Lower New York River
           Inst != "CHOPAX_ChespkBay_Jellyfish_Data" & # Chesapeake Bay
           Inst != "TIES_ChespkBay" & # Chesapeake Bay
           Inst != "ChespkBay_Prog" & # Chesapeake Bay
           Inst != "PortPhillipBay" & # Port Phillip Bay
           Inst != "NarragansettBay_Plankton" & # In Rhode Island Sound
           Inst != "Norwegian_Literature" & 
           Inst != "Dutch_Literature" & # Chesapeake Bay
           Inst != "Salps_NSW_AUS" &
           Inst != "Salp_Tasman_NSW_AUS") %>% 
  mutate_at(vars(rank_phylum), list(~ifelse(taxon == "medusa", "Cnidaria", .)))

jedi <- jedi %>%
  mutate(Group = rank_phylum,
         Group = str_replace(Group,"Cnidaria","Jellyfish"),
         Group = str_replace(Group,"salp","Salps"),
         Group = as.factor(Group),
         TimeLocal = hour(TimeLocal) + (minute(TimeLocal) + (second(TimeLocal)/60))/60,
         IDforSumming = paste0(sprintf('%02.2f', abs(Latitude)),
                               sprintf('%03.2f', abs(Longitude)),
                               sprintf('%4d', Year),
                               sprintf('%2d', Month),
                               sprintf('%2d', Day)), # Create Unique ID for every net deployment
         IDforSumming = as.factor(IDforSumming)) %>% 
  droplevels()


# Now I have the data, I need to fix up missing data, mesh, time etc
jedi <- jedi %>%
  mutate_at(vars(Mesh), list(~ifelse(Inst == "Helgoland_Roads", 500, .))) %>% # Helgolands uses 2 nets, but I assume this data is the macrozoop (500um) rather than mesozoop (hence 150um)
  mutate_at(vars(Upper_Z), list(~ifelse(Inst == "Helgoland_Roads", 1, .))) %>% 
  mutate_at(vars(Depth), list(~ifelse(Inst == "Helgoland_Roads", 8, .))) %>% 
  mutate_at(vars(Lower_Z), list(~ifelse(Inst == "Helgoland_Roads", 4, .))) %>% 
  mutate_at(vars(Tow), list(~ifelse(Inst == "Helgoland_Roads", "O", .))) %>% 
  mutate_at(vars(TimeLocal), list(~ifelse(Inst == "Helgoland_Roads", runif(length(Inst), min = 6, max = 9), .))) # Samples were generally taken before 9 am Wiltshire, K.H. & Dürselen, CD. Helgol Mar Res (2004) 58: 252. https://doi.org/10.1007/s10152-004-0192-4

jedi <- jedi %>%
  mutate_at(vars(Mesh), list(~ifelse(Inst == "PACES", 76, .))) %>% 
  mutate_at(vars(Depth), list(~ifelse(Inst == "PACES", 0, .))) %>% 
  mutate_at(vars(Upper_Z), list(~ifelse(Inst == "PACES", 0, .))) %>% 
  mutate_at(vars(Lower_Z), list(~ifelse(Inst == "PACES", 0, .)))

# Citation: Escribano, R., Hidalgo, P., Fuentes, M. and Donoso, K., 2012. Zooplankton 
# time series in the coastal zone off Chile: variation in upwelling and responses of 
# the copepod community. Progress in Oceanography, 97, pp.174-186.
# This time series appears to be the St 18 Concepcion.
jedi <- jedi %>%
  mutate_at(vars(Mesh), list(~ifelse(Inst == "HumboldtCurrent_Survey", 200, .))) %>% # Helgolands uses 2 nets, but I assume this data is the macrozoop (500um) rather than mesozoop (hence 150um)
  mutate_at(vars(Upper_Z), list(~ifelse(Inst == "HumboldtCurrent_Survey", 1, .))) %>% 
  mutate_at(vars(Depth), list(~ifelse(Inst == "HumboldtCurrent_Survey", 40, .))) %>% 
  mutate_at(vars(Lower_Z), list(~ifelse(Inst == "HumboldtCurrent_Survey", 80, .))) %>% 
  mutate_at(vars(Tow), list(~ifelse(Inst == "HumboldtCurrent_Survey", "O", .))) %>% 
  mutate_at(vars(Gear), list(~ifelse(Inst == "HumboldtCurrent_Survey", 176, .))) %>% # COPEPOD Code for Tucker Trawl
  mutate_at(vars(TimeLocal), list(~ifelse(Inst == "HumboldtCurrent_Survey", runif(length(Inst), min = 12, max = 16), .))) # Samples were generally taken before 9 am Wiltshire, K.H. & Dürselen, CD. Helgol Mar Res (2004) 58: 252. https://doi.org/10.1007/s10152-004-0192-4



# Citation: Dutto MS, Chazarreta CJ, Rodriguez CS, Schiariti A, Diaz Briz LM, Genzano GN 
# (2019) Macroscale abundance patterns of hydromedusae in the temperate Southwestern 
# Atlantic (27 ̊–56 ̊ S). PLoS ONE 14(6): e0217628. https://doi.org/ 10.1371/journal.pone.0217628

jedi <- jedi %>%
  mutate_at(vars(Mesh), list(~ifelse(Inst == "Proyecto_Anchoita_INIDEP", 350, .))) %>% 
  mutate_at(vars(TimeLocal), list(~ifelse(Inst == "Proyecto_Anchoita_INIDEP", runif(length(Inst), min = 6, max = 18), .))) %>% 
  mutate(Group = str_replace(Group, "Chordata","Salps"),
         Group = as.factor(Group))


### I NEED TO DO THIS IN ORDER TO GET THE SAMPLE DEPTHS FOR INIDEP WHICH WENT FROM THE BOTTOM.
########### IMPORT BATHYMETRY ###########
EnviroDir <- "Data/EnvironmentalData"
bathy_files <- "GEBCO_BATHY_2002-01-01_rgb_360x180.csv"
bathy_raw <- read.csv(paste0(EnviroDir,"/",bathy_files), header = FALSE)
bathy_vector <- as.vector(t(as.matrix(bathy_raw))) # Convert to vector

# create a grid of lat/lon, add headers, and convert to data frame
bathy_matrix = matrix(NA, nrow = length(bathy_vector), ncol = 3)
bathy_matrix[,c(1,2)] = as.matrix(expand.grid("Long" = seq(-179.5, 179.5, 1), "Lat" = seq(89.5, -89.5, -1)))
bathy_matrix[,3] <- bathy_vector
colnames(bathy_matrix) <- c("Long", "Lat", "Bathy")
bathy_matrix <- as.data.frame(bathy_matrix)
bathy_matrix[c(bathy_matrix$Bathy > 0), "Bathy"] <- NA
bathy_matrix <- bathy_matrix[!is.na(bathy_matrix$Bathy),] 
nvb <- ann(as.matrix(bathy_matrix[,c("Long", "Lat")]), as.matrix(jedi[,c("Longitude", "Latitude")]), 
           k = 1, verbose = FALSE)$knnIndexDist[,1]
jedi$Bathy <- bathy_matrix[nvb,3]

rm(bathy_files, bathy_raw, bathy_vector, bathy_matrix, EnviroDir, nvb)

jedi <- jedi %>%
  mutate(Bathy = abs(Bathy),
         Bathy = case_when(Bathy < 10 ~ 10,
                           Bathy >= 10 ~ Bathy)) %>% 
  mutate_at(vars(Tow), list(~ifelse(Inst == "Proyecto_Anchoita_INIDEP", "O", .))) %>% 
  mutate_at(vars(Upper_Z), list(~ifelse(Inst == "Proyecto_Anchoita_INIDEP", 0, .))) %>% 
  mutate_at(vars(Lower_Z), list(~ifelse(Inst == "Proyecto_Anchoita_INIDEP", Bathy, .))) %>% 
  mutate_at(vars(Depth), list(~ifelse(Inst == "Proyecto_Anchoita_INIDEP", (Lower_Z-Upper_Z)/2, .))) %>% 
  select(-Bathy) # Now remove Bathy bectause I don't want it here yet.

  
  # Calculate total abundance for each sample (IDforSumming)
jedi2 <- jedi %>% group_by(Group, IDforSumming) %>%
  summarise(TotAbundance = sum(Abundance))

jedi <- jedi %>% 
  ungroup(jedi) %>% 
  droplevels()

# All rows in a IDforSumming are the same except for ScientificName, so keep the distinct ones
jedi <- jedi %>% mutate(Group = as.factor(Group)) %>%
  group_by(Group, IDforSumming) %>%
  distinct(Group, IDforSumming, .keep_all = TRUE) %>% 
  right_join(jedi2, by = c("Group", "IDforSumming")) %>% 
  ungroup() %>% 
  dplyr::select(Retain) %>% 
  droplevels()

rm(jedi2) # Clean up

hist(jedi$TimeLocal) # Problem with the JeDI time distribution. I suspect many day or night samples have been placed at 12 am/pm
summary(jedi)

saveRDS(jedi,'JeDI_Final.rds')
