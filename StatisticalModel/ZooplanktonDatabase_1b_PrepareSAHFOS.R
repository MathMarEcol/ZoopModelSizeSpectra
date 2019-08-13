rm(list = ls())

library(tidyverse)
library(lubridate)
library(worrms)
re_dld <-  0

Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")

######################### Read and clean up data ####################################
sahfos <- read_csv("Data/SAHFOS/raw_cpr_survey_extract.csv")

sahfos <- sahfos %>% 
  rename(Latitude = sample_latitude, Longitude = sample_longitude,
         Abundance = count_per_m3, aphiaID = aphia_id, 
         UniqueSampleID = Sample_ID) %>%
  mutate(TimeGMT = hour(UTC_Date) + (minute(UTC_Date) + (second(UTC_Date)/60))/60,
         Year = year(UTC_Date),
         Month = month(UTC_Date),
         Day = day(UTC_Date),
         TimeLocal = TimeGMT + round(Longitude/15),
         ShipCruise = str_extract(UniqueSampleID,"^[^\\-]*"),
         Transect = str_extract(ShipCruise,"\\D+"),
         Transect = as.factor(Transect)) %>% 
  add_column(Type = "CPR", DBase = "SAHFOS", Depth = NA, Tow = "S",
         Upper_Z = NA, Lower_Z = NA, Project = "SAHFOS", Gear = 191, 
         Mesh = 270, Volume = NA, Ship = NA) %>% 
  mutate(TimeLocal = case_when(TimeLocal >= 0 & TimeLocal <= 24 ~ TimeLocal,
                               TimeLocal > 24 ~ TimeLocal - 24,
                               TimeLocal < 0 ~ TimeLocal + 24)) %>% 
  select(-UTC_Date) %>%
  filter(!str_detect(taxon_name, "Calanus Total Traverse") &
           !str_detect(taxon_name, "Calanus TOTAL \\(total count\\)") &
           !str_detect(taxon_name, "Centropages chierchiae traverse") &
           !str_detect(taxon_name, "Copepod eggs") &
           !str_detect(taxon_name, "Copepoda TOTAL \\(total count\\) \\-\\[ \\]\\-") &
           !str_detect(taxon_name, "Metridia Total traverse") &
           !str_detect(taxon_name, "Euphausiacea Total") &
           # !str_detect(taxon_name, "Euphausiacea \\-\\[ eggs \\]\\-") &
           !str_detect(taxon_name, "Chaetognatha Total traverse") &
           !str_detect(taxon_name, "Thaliacea") &
           !str_detect(taxon_name, "Plastic") &
           !str_detect(taxon_name, "Parasites of the plankton") &  
           !str_detect(taxon_name, "Fish eggs") &  
           !str_detect(taxon_name, "Stellate body") &  
           !str_detect(taxon_name, "Spiny egg") &  
           !str_detect(taxon_name, "Clione shells")) %>% 
  mutate(aphiaID = as.numeric(aphiaID),
         taxon_name = str_trim(taxon_name),
         taxon_name = str_remove(taxon_name, " \\(\\'Exuviaella\\' type\\)"),
         taxon_name = str_remove(taxon_name, " \\(Antarctica\\)"),
         taxon_name = str_remove(taxon_name, " \\(Unidentified\\)"),
         taxon_name = str_remove(taxon_name, " Adult"),
         taxon_name = str_remove(taxon_name, " Dinoflagellate cysts"),
         taxon_name = str_remove(taxon_name, " eyecount"),
         taxon_name = str_remove(taxon_name, " eyecount"),
         taxon_name = str_remove(taxon_name, " I-V"),
         taxon_name = str_remove(taxon_name, " Total Traverse"),
         taxon_name = str_remove(taxon_name, " total traverse"),
         taxon_name = str_remove(taxon_name, " Total"),
         taxon_name = str_remove(taxon_name, " VI Female"),
         taxon_name = str_remove(taxon_name, " VI Female"),
         taxon_name = str_remove(taxon_name, " VI Male"),
         taxon_name = str_remove(taxon_name, regex("(\\s\\&).*")),
         taxon_name = str_remove(taxon_name, regex("(\\s\\().*")),
         taxon_name = str_remove(taxon_name, regex("(\\sV-VI).*")),
         taxon_name = str_remove(taxon_name, regex("(\\VI).*")),
         taxon_name = str_remove(taxon_name, regex("(\\sII).*")),
         taxon_name = str_remove(taxon_name, regex("(\\sCVI).*")),
         taxon_name = str_remove(taxon_name, regex("(\\sV).*")),
         taxon_name = str_remove(taxon_name, regex("(\\sC).*")),
         taxon_name = str_remove(taxon_name, regex("\\..*")),
         taxon_name = str_remove(taxon_name, regex("\\(.*")),
         taxon_name = str_remove(taxon_name, regex("\\(\\sI\\-V\\)\\.*")),
         taxon_name = str_remove(taxon_name, regex("\\(\\sV\\-VI\\)\\.*")),
         taxon_name = str_replace(taxon_name, " \"CPR-TOTAL \"\"eyecount\"\"\"",""),
         taxon_name = str_replace(taxon_name, " \"CPR-TOTAL \"\"traverse\"\"\"",""),
         taxon_name = str_replace(taxon_name, " parts/fragments",""),
         taxon_name = str_replace(taxon_name, " sp\\.", ""),
         taxon_name = str_replace(taxon_name, " sp\\.", ""),
         taxon_name = str_replace(taxon_name, " spp", ""),
         taxon_name = str_replace(taxon_name, " spp\\.", ""),
         taxon_name = str_replace(taxon_name, " spp\\.", ""),
         taxon_name = str_replace(taxon_name, " TOTAL \\(total count\\)",""),
         taxon_name = str_replace(taxon_name, "Copepod Eggs", "Copepoda"),
         taxon_name = str_replace(taxon_name,"Small Calanoid copepod \\(Unidentified\\)","Calanus"),
         taxon_name = str_replace(taxon_name,"Para Pseudocalanus","Pseudocalanus")
  )

######################### Get aphiaID and classify all groups ####################################
lookup <- sahfos %>% 
  select(c(taxon_name, aphiaID)) %>%
  distinct(aphiaID, .keep_all = TRUE) %>%
  mutate(aphiaID = as.numeric(aphiaID))

## Get the classification information for the SAHFOS data in order to attribute the classifcation. 
## We already have the aphiaID for this dataset - Just need Copepod etc
if (re_dld == 1){
  out <- wm_classification_(id = lookup$aphiaID)
  out$id <- as.numeric(out$id)
  saveRDS(out,"Data/SAHFOS/Aphia_SAHFOS.rds")
}else{
  out <- readRDS("Data/SAHFOS/Aphia_SAHFOS.rds")
}

# Now create a lookup table - This is really clunky. There must be a better way.....but it works
ID <- out$id[out$scientificname=="Copepoda"] # These are the species' IDs
aphiaID <- data.frame("aphiaID"=ID,stringsAsFactors=FALSE)
aphiaID$Group <- "Copepods"
rm(ID)

ID <- out$id[out$scientificname=="Chaetognatha"]
temp <- data.frame("aphiaID"=ID,stringsAsFactors=FALSE)
temp$Group <- "Chaetognaths"
aphiaID <- bind_rows(aphiaID,temp)
rm(temp,ID)

ID <- out$id[out$scientificname=="Euphausiacea"]
temp <- data.frame("aphiaID"=ID,stringsAsFactors=FALSE)
temp$Group <- "Euphausiids"
aphiaID <- bind_rows(aphiaID,temp)
rm(temp,ID)

ID <- out$id[out$scientificname=="Scyphozoa" | out$scientificname=="Hydrozoa"]
temp <- data.frame("aphiaID"=ID,stringsAsFactors=FALSE)
temp$Group <- "Jellyfish"
aphiaID <- bind_rows(aphiaID,temp)
rm(temp,ID)

ID <- out$id[out$scientificname=="Salpida"]
temp <- data.frame("aphiaID"=ID,stringsAsFactors=FALSE)
temp$Group <- "Salps"
aphiaID <- bind_rows(aphiaID,temp)
rm(temp,ID)

# Appendicularia - Need to check for both
# wm_classification_(id = 17446) # old
# wm_classification_(id = 146421) # new
ID <- out$id[out$scientificname=="Appendicularia" | out$scientificname=="Larvacea"]
temp <- data.frame("aphiaID"=ID,stringsAsFactors=FALSE)
temp$Group <- "Larvaceans"
aphiaID <- rbind(aphiaID,temp)

rm(lookup)

# Create a unique list of names which are unmatched before I move on.
sahfos_na <- sahfos %>% 
  left_join(aphiaID,by="aphiaID") %>%
  filter(is.na(Group)) %>% 
  distinct(taxon_name) %>% 
  arrange(taxon_name)

sahfos <- sahfos %>% 
  left_join(aphiaID,by="aphiaID") %>%
  filter(!is.na(Group)) %>%
  mutate(aphiaID = as.numeric(aphiaID)) %>% 
  select(-taxon_name)

#############################################
# Check species names the same in sahfos$SpeciesName and FeedingType$TAXON_NAME - Yes!!!
# di <- dplyr::setdiff(levels(FeedingType$TAXON_NAME), levels(sahfos$taxon_name))

FeedingType <- read_csv("Data/Copepod_Feeding_Type_v3_wAPHIA.csv")
FeedingType <- FeedingType %>% 
  select(-c(TAXON_NAME, GENUS_NAME)) %>%
  distinct()

## This next line has increased the size of sahfos. 
# I need to work out why and sdtop it. It is 
# replicating lines for some reason. It is joining because 
# there are no common IDs between the tibbles. Most of the 
# sahfos alphaIDs are species, whereas all the feeding ones are genus

# These are the GENUS AphiaIDs. Now to merge with SAHFOS so we can merge with Feeding
# Unfortuantely this implementation won't select the Family level matches
# But probably fine for the moment.
Genus <- out %>% filter(rank == "Genus") %>% # AphiaID is the ID of the GENUS
  select(-c(rank,scientificname)) %>%
  rename(GenusID = AphiaID, aphiaID = id)

sahfos2 <- sahfos %>% 
  left_join(Genus,GenusID, by = "aphiaID") %>% # Add the Genus ID so I can merge with FeedingType
  left_join(FeedingType, by = c("GenusID" = "aphiaID")) %>% # Merge with feeding type
# Set Diets: All copepod species assigned a diet
  mutate(Diet = case_when(FEED == "CC" ~ "Carnivore", # Diet: Carnivores (CC) vs 
                          FEED == "CBF" ~ "Omnivore", # Omnivores (CBF + CD + CH + CO + CP + SF)
                          FEED == "CD" ~ "Omnivore",
                          FEED == "CH" ~ "Omnivore",
                          FEED == "CO" ~ "Omnivore",
                          FEED == "CP" ~ "Omnivore",
                          FEED == "SF" ~ "Omnivore"),
         Diet = as.factor(Diet),
         Group = replace(Group, Diet=="Carnivore", "CarnCopepods"),
         Group = replace(Group, Diet=="Omnivore", "OmniCopepods")) %>%
  mutate(Diet2 = case_when(FEED == "CC" ~ "Carnivore", # Diet2: vs Carnivore (CC + CBF + CD + CH + CO + CP + SF)
                           FEED == "CBF" ~ "Carnivore",
                           FEED == "CD" ~ "Carnivore",
                           FEED == "CH" ~ "Carnivore",
                           FEED == "CO" ~ "Carnivore",
                           FEED == "CP" ~ "Carnivore",
                           FEED == "SF" ~ "Carnivore",
                           FEED == "CH" ~ "Herbivore"), # Herbivore (CH) 
         Diet2 = as.factor(Diet2)) %>%
  filter(Group!="Copepods")

# # Not all spcies are matched - Need to fix this
# un_matched <- sahfos2 %>% filter(Group=="Copepods") %>%
#   distinct(taxon_name)

####################################################

sahfos_wide <- sahfos2 %>% 
  spread(key = Group, value = Abundance, fill = 0) %>%  # Wide format: Move the Abundance into each group.
  select(-c("GenusID", "FEED", "Diet", "Diet2", "Taxon_ID","aphiaID")) %>% 
  arrange(UniqueSampleID, Year, Month, Day) %>%  # Arrange to make it easier to debug
  group_by(UniqueSampleID) # Regroup by UniqueSampleID


# Calculate total abundance for each UniqueSampleID
sahfos_wide2 <- sahfos_wide %>% 
  summarise(Chaetognaths = sum(Chaetognaths),
            CarnCopepods = sum(CarnCopepods),
            OmniCopepods = sum(OmniCopepods),
            Euphausiids = sum(Euphausiids),
            Jellyfish = sum(Jellyfish),
            Larvaceans = sum(Larvaceans),
            Salps = sum(Salps))

# Remove the Groups from sahfos_wide so we don't get duplicates on the inner_join
sahfos_wide3 <- sahfos_wide %>% select(-c("Chaetognaths", "CarnCopepods", "OmniCopepods", 
                                          "Euphausiids", "Jellyfish", "Larvaceans", "Salps")) %>%
  inner_join(sahfos_wide2, by="UniqueSampleID") # Join sahfos_wide (metadata) with sahfos_wide2 (Abundance)

sahfos_final <-  sahfos_wide3 %>% 
  gather(c("Chaetognaths", "CarnCopepods", "OmniCopepods", "Euphausiids", "Jellyfish", "Larvaceans", "Salps"), 
         key= "Group",value = "TotAbundance") %>% 
  arrange(UniqueSampleID, Year, Month, Day, TimeGMT) %>%  # Arrange to make it easier to debug
  distinct(Group, UniqueSampleID, .keep_all = TRUE) %>% 
  select(Retain) %>% 
  ungroup()

saveRDS(sahfos_final,'DatabaseOutput/SAHFOS_Final.rds')

rm(sahfos, sahfos_na, sahfos2, temp, sahfos_wide, sahfos_wide2, sahfos_wide3, aphiaID, FeedingType, Genus, out)
