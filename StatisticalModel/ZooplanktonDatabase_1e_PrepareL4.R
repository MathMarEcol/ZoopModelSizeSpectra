rm(list = ls())

library(tidyverse)
library(lubridate)
# Station L4 (50 15 N, 48 13 W) is located waters 7.5 nautical miles (􏰁13.9km)
# south-west off Plymouth. Since March 1988, zooplankton samples have been collected 
# weekly, weather permitting, by ver- tical net hauls from the sea floor to the surface 
# using a WP2 net with a mesh-size of 200 mm and a mouth area of 0.25 m2 (UNESCO, 1968). 
# Two hauls are succes- sively taken at approximately mid-morning and the samples are 
# preserved and stored in 5% formalin. Zooplankton are later sub-sampled, counted and 
# ident- ified under a microscope in the laboratory. Subsamples are extracted using a 
# Folsom splitter and a Stempel pipette, to identify separately large and small organisms. 
# Subsamples contain around 200 – 400 individuals. Abundances in the two hauls taken on 
# each sampling date are averaged to reduce the variability related to the sampling. 
# Abundance is expressed as numbers of organ- isms per cubic meter (N m23).


Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")

L4 <- read_csv("Data/Plymouth_L4_1988-2017_ZoopAbund_v15.2.19.csv")

L4 <- L4 %>% 
  add_column(Volume=NA, Depth=NA, Upper_Z=0, Lower_Z=55, 
             Tow="V", Mesh=200, Ship=NA, Group=NA, TimeLocal=10, TimeGMT = NA,
             ShipCruise="L4", DBase="L4", Gear=132, Type="Net", Project = "L4",
             Transect="L4") %>% # Add missing variables
  rename(Date = SamplingDate, UniqueSampleID = DateOrder) %>% 
  mutate(Date = dmy(Date),
         Year = year(Date),
         Month = month(Date),
         Day = day(Date),
         ShipCruise = as.factor(ShipCruise),
         Day = as.numeric(Day),
         Volume = as.numeric(Volume),
         Transect = as.factor(Transect)) %>% 
  dplyr::select(-c(Season))

L4_long <- L4 %>% gather(Taxa,Abundance,"Anemone larvae":"Copepod nauplii")

Feed <- read_csv("Data/L4_Feeding.csv")
L4_long <- L4_long %>% left_join(Feed, by="Taxa")

L4 <- L4_long %>% 
  mutate(Group = replace(Group, Feed=="Omnivore",'OmniCopepods'),
         Group = replace(Group, Feed=="Carnivore",'CarnCopepods'),
         Group = replace(Group, Taxa=="Total medusae",'Jellyfish'),
         Group = replace(Group, Taxa=="Total Chaetognath",'Chaetognaths'),
         Group = replace(Group, Taxa=="Euphausiid nauplii",'Euphausiids'),
         Group = replace(Group, Taxa=="Euphausiid calyptopis",'Euphausiids'),
         Group = replace(Group, Taxa=="Euphausiid furcilia",'Euphausiids'),
         Group = replace(Group, Taxa=="Euphausiid adult",'Euphausiids'),
         Group = replace(Group, Taxa=="Appendicularia",'Larvaceans')) %>% 
  filter(!is.na(Group)) %>% 
  group_by(Group, UniqueSampleID)

########  Create totals ###############
L4sum <- L4 %>% 
  dplyr::summarise(TotAbundance = sum(as.numeric(Abundance)))# Calculate total abundance for each UniqueSampleID

# All rows in a UniqueSampleID are the same except for ScientificName, so keep the distinct ones 
L4 <- L4 %>% distinct(Group, UniqueSampleID, .keep_all = TRUE) %>% 
  right_join(L4sum, by = c("Group", "UniqueSampleID")) %>%
  ungroup() %>% 
  mutate(Group = as.factor(Group),
         UniqueSampleID = as.factor(UniqueSampleID)) %>% 
  dplyr::select(Retain)
                    

L4wide <- spread(L4, Group, TotAbundance, fill = 0) # Wide format: Move the Abundance into each group. 
L4wide <- arrange(L4wide, ShipCruise, Year, Month, Day, TimeGMT) # Arrange to make it easier to debug
L4wide <- group_by(L4wide, UniqueSampleID) # Regroup by UniqueSampleID

# Calculate total abundance for each UniqueSampleID
L4wide2 <- summarise(L4wide,
                       Chaetognaths = sum(Chaetognaths),
                       CarnCopepods = sum(CarnCopepods),
                       OmniCopepods = sum(OmniCopepods),
                       Euphausiids = sum(Euphausiids),
                       Jellyfish = sum(Jellyfish),
                       Larvaceans = sum(Larvaceans))
# Remove the Groups from CPR_wide so we don't get duplicates on the inner_join
L4wide <- dplyr::select(L4wide,-c("Chaetognaths", "CarnCopepods", "OmniCopepods", "Euphausiids", "Jellyfish", "Larvaceans"))
L4wide3 <- inner_join(L4wide, L4wide2, by="UniqueSampleID") # Join CPR_wide (metadata) with CPR_wide2 (Abbundance)

L4long <- gather(L4wide3, Group, TotAbundance, Chaetognaths:Larvaceans, factor_key = TRUE)
L4long <- arrange(L4long, ShipCruise, Year, Month, Day, TimeGMT) # Arrange to make it easier to debug

L4long <- ungroup(L4long)


# save data
saveRDS(L4long,'L4_Final.rds')


