library(tidyverse)
library(taxize)
library(worrms)

# It is unclear what the difference is between TAXON_NAME and NEW_NAME
re_dld <- 1 # Should we download the data

###### Import in FeedingType and clean it and match to aphiaID
FeedingType <- read_csv("Data/Copepod_Feeding_Type_v3.csv")

# Cyclops Feeding to CP because there are two version in the Feeding Spreadsheet
# CHECK WITH ANT/CLAIRE EVENTUALLY
FeedingType$FEED[FeedingType$GENUS_NAME=="Cyclops"] = "CP"

# Reduce to the min numer of unique names
FeedingType_reduce <- FeedingType %>%  distinct(GENUS_NAME, .keep_all = TRUE) %>%
  select(c(GENUS_NAME))

# Get classification ID
if (re_dld==1){
  # Now extract the highest classification and the ID so we can match back with FeedingType
  FeedingType_reduce$name2 <- NA # Name2 is just for checking the match up
  FeedingType_reduce$aphiaID <- NA
  for (i in 1:length(FeedingType_reduce$GENUS_NAME)){
    aphiaID <- wm_records_taxamatch(FeedingType_reduce$GENUS_NAME[i], marine_only = FALSE)
    FeedingType_reduce$name2[i] <- aphiaID[[1]]$scientificname
    FeedingType_reduce$aphiaID[i] <- aphiaID[[1]]$AphiaID
    rm(aphiaID)
  }
  saveRDS(FeedingType_reduce, file = "Data/aphiaID_4FeedingTypeGENUS.rds")
} else{
  FeedingType_reduce = readRDS("Data/aphiaID_4FeedingTypeGENUS.rds")
}

########################################

# Then we joing to FeedingType
FeedingType <- FeedingType %>% 
  left_join(FeedingType_reduce,by="GENUS_NAME")
rm(FeedingType_reduce)

# Check the match up works, then delete the name2 column
FeedingType <- FeedingType %>% select(-c(name2))

# Size is 0 - great. No unamtched Genus
FT_NA <- dplyr::filter(FeedingType,is.na(aphiaID)) 

# But there are unmatched Feeding Types (from the input data). 
# Remove these
FeedingType <- FeedingType %>% 
  filter(!is.na(FEED)) 

write_csv(FeedingType,"Data/Copepod_Feeding_Type_v3_wAPHIA.csv")




########## OLD CODE - NOT AS MUCH CLEANING NEEDED NOW I USE GENUS ############
# # Clean up species names to allow automatic retrieval of aphiaID
# FeedingType <- FeedingType %>% mutate(aphiaName = GENUS_NAME,
#                                       aphiaName = str_replace(aphiaName, " spp\\.", ""),
#                                       aphiaName = str_replace(aphiaName, " sp\\.", ""),
#                                       aphiaName = str_replace(aphiaName, " \"CPR-TOTAL \"\"traverse\"\"\"",""),
#                                       aphiaName = str_replace(aphiaName, " TOTAL \\(total count\\)",""),
#                                       aphiaName = str_replace(aphiaName, " \"CPR-TOTAL \"\"eyecount\"\"\"",""),
#                                       aphiaName = str_remove(aphiaName, regex("(\\s\\&).*")),
#                                       aphiaName = str_replace(aphiaName, " parts/fragments",""))

# There are errors in the Feeding Database. I will let Claire/Ant sort out what is wrong, 
# in case these changes need to go right through the CSIRO database
# FeedingType <- FeedingType %>% mutate(aphiaName = str_replace(aphiaName, "Arietellus armata", "Arietellus armatus"),
#                                       aphiaName = str_replace(aphiaName, "Centropages elongata", "Centropages elongatus"),
#                                       aphiaName = str_replace(aphiaName, "Centropages kroeyeri", "Centropages kroyeri"),
#                                       aphiaName = str_replace(aphiaName, "Centropages orsini", "Centropages orsinii"),
#                                       aphiaName = str_replace(aphiaName, "Chiridiella abussalis", "Chiridiella abyssalis"),
#                                       aphiaName = str_replace(aphiaName, "Clausocalanus breuipes", "Clausocalanus brevipes"),
#                                       aphiaName = str_replace(aphiaName, "Clausocalanus major", "Clausocalanus arcuicornis major"),
#                                       aphiaName = str_replace(aphiaName, "Cyclopina litoralis", "Cyclopinoides littoralis"),
#                                       aphiaName = str_replace(aphiaName, "Diaixis pigmaea", "Diaixis pygmaea"),
#                                       aphiaName = str_replace(aphiaName, "Euaugaptilus similis", "Euaugaptilus similis"),
#                                       aphiaName = str_replace(aphiaName, "Eucalanus pseudoattenuatus", "Eucalanus pseudattenuatus"),
#                                       aphiaName = str_replace(aphiaName, "Euchirella mesinensis", "Euchirella messinensis"),
#                                       aphiaName = str_replace(aphiaName, "Gaetanus brevicaudatus", "Gaetanus brevicaudatus"),
#                                       aphiaName = str_replace(aphiaName, "Macrosetella sulcata", "Macrosetella"),
#                                       aphiaName = str_replace(aphiaName, "Mecynocera Microcalanus", ""),
#                                       aphiaName = str_replace(aphiaName, "Pareuchaeta rasa", "Paraeuchaeta rasa"),
#                                       aphiaName = str_replace(aphiaName, "Pleuromamma robusta f antarctica", "Pleuromamma robusta antarctica"),
#                                       aphiaName = str_replace(aphiaName, "Scaphocalanus minutus", "Scaphocalanus minuta"),
#                                       aphiaName = str_replace(aphiaName, "Scolecithricella gracialis", "Scolecithricella glacialis"),
#                                       aphiaName = str_replace(aphiaName, "Scolecithrichopsis ctenopus", "Scolecitrichopsis ctenopus"))

