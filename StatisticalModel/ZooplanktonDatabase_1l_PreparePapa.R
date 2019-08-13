# Waddell, Brenda J., and Skip McKinnell. 1995. Ocean Station "Papa" 
# detailed zooplankton data: 1956- 1980. Can. Tech. Rep. Fish. Aquat. 
# Sci. 2056: 21 p.

# Zooplankton samples were collected at Ocean Station "P" (50°N, 145°W) 
# from 1956 to 1980, and were analyzed to various levels of taxonomic 
# resolution over the years. Although summaries of these data have been 
# previously published by LeBrasseur (1965) and Fulton (1978, 1983), 
# the detailed species data have never been published. We have reformatted 
# the detailed data, corrected any errors we discovered, and added extra 
# information to produce one complete dataset for the whole sampling period. 
# This dataset contains total zooplankton wet weights/m3 for the whole 
# period 1956 to 1980, as well as densities (numbers/m3) for five major 
# taxa (copepods, chaetognaths, euphausiids, amphipods, and Aglantha) from 
# 1964 to 1967, and species identifications, counts and lengths for many 
# samples collected between 1968 to 1980. The purpose of this document is 
# to make the detailed data available to the scientific community in 
# electronic format, and to provide a convenient reference for citing the 
# detailed data.

library(tidyverse)
library(lubridate)

Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")

Papa <- read_csv("Data/OceanStationPapa/Zooplankton_1956-1981_OSPZOO_EN.csv",
                 col_types = "dccdDcddcddcddd", na = c("", "."))

Papa <- Papa %>% 
  drop_na(Density) %>% 
  rename(ShipCruise = `Cruise Number`, SampleID = `Sample Number`, TimeLocal = `Time of Day`,
         Date = `Date (yyyy-mm-dd)`, Lower_Z = `Sample Depth (m)`, Abundance = Density, 
         CODE = `Species Code (alpha)`) %>% 
  select(-c(`Gear (alpha)`, `Technician Initials (alpha)`, `Wire Angle (o)`, `Sample Analysis Status (alpha)`,
            `Sample Weight (g)`, `Wet Weight Biomass (mg/m3)`, `Size (mm)`, `Count`))

# Fix up local time
Papa <- Papa %>% 
  mutate(TimeLocal = as.numeric(TimeLocal),
         H = floor(TimeLocal/100),
         M = TimeLocal-H*100,
         TimeLocal = H + M/60) %>% 
  select(-c(H, M)) %>% 
  mutate(Year = year(Date),
         Month = month(Date),
         Day = day(Date)) %>% 
  select(-Date)

species <- read_csv("Data/OceanStationPapa/OSP_SpeciesCodes.csv", na = c("", "."), skip = 4)

species <- species %>% 
  mutate(CODE = str_replace(CODE, "RKO","RK0")) # Replace O with zero (0)

Papa <- left_join(Papa, species, by = "CODE")

Papa <- Papa %>% 
  drop_na(SPECIES) %>% 
  mutate(Genus = str_extract(SPECIES,"^\\S*"),
         Genus = str_replace(Genus, "Harpacticoid","Harpacticoida"),
         Genus = str_replace(Genus, "Harpacticoidaa","Harpacticoida"),
         Genus = str_replace(Genus, "Oikpleura","Oikopleura"))

## Load and check the Copepod Feeding
FeedingType <- read_csv("Data/Copepod_Feeding_Type_v3_wAPHIA.csv")
FeedingType <- FeedingType %>% 
  dplyr::select(-c(TAXON_NAME, aphiaID)) %>%
  distinct()

Papa <- Papa %>% 
  left_join(FeedingType, by = c("Genus" = "GENUS_NAME")) %>% 
  mutate(Group = case_when(FEED == "CC" ~ "CarnCopepods", # Diet: Carnivores (CC) vs 
                           FEED == "CBF" ~ "OmniCopepods", # Omnivores (CBF + CD + CH + CO + CP + SF)
                           FEED == "CD" ~ "OmniCopepods",
                           FEED == "CH" ~ "OmniCopepods",
                           FEED == "CO" ~ "OmniCopepods",
                           FEED == "CP" ~ "OmniCopepods",
                           FEED == "SF" ~ "OmniCopepods",
                           str_detect(Genus,"Euphausiacea") ~ "Euphausiids",
                           str_detect(Genus,"Euphausia") ~ "Euphausiids",
                           str_detect(Genus,"Nematoscelis") ~ "Euphausiids",
                           str_detect(Genus,"Thysanoessa") ~ "Euphausiids",
                           str_detect(Genus,"Tessarabrachion") ~ "Euphausiids",
                           str_detect(Genus,"Nematoscelis") ~ "Euphausiids",
                           str_detect(Genus,"Salpa") ~ "Salps",
                           str_detect(Genus,"Thalia") ~ "Salps",
                           str_detect(Genus,"Chaetognatha") ~ "Chaetognaths",
                           str_detect(Genus,"Sagitta") ~ "Chaetognaths",
                           str_detect(Genus,"Eukrohnia") ~ "Chaetognaths",
                           str_detect(Genus,"Larvacea") ~ "Larvaceans",
                           str_detect(Genus,"Oikopleura") ~ "Larvaceans",
                           str_detect(Genus,"Fritillaria") ~ "Larvaceans",
                           str_detect(Genus,"Medusa") ~ "Jellyfish",
                           str_detect(Genus,"Aglantha") ~ "Jellyfish",
                           str_detect(Genus,"Proboscidactyla") ~ "Jellyfish",
                           str_detect(Genus,"Rathkea") ~ "Jellyfish",
                           str_detect(Genus,"Hybocodon") ~ "Jellyfish",
                           str_detect(Genus,"Velella") ~ "Jellyfish",
                           str_detect(Genus,"Agalma") ~ "Jellyfish",
                           str_detect(Genus,"Dimophyes") ~ "Jellyfish",
                           str_detect(Genus,"Chelophyes") ~ "Jellyfish",
                           str_detect(Genus,"Muggiaea") ~ "Jellyfish"),
                           Group = as.factor(Group)) %>% 
           filter(!is.na(Group)) %>% 
           select(-c(SPECIES, Genus, URI, CODE, FEED))


Papa <- Papa %>%      
  add_column(Volume=NA, Upper_Z=0, Tow="V", Mesh=351, DBase="OSP", 
             Type="Net", Project = "OSP", Transect="OSP", 
             Latitude = 50, Longitude = -145, TimeGMT = NA) %>% # Add missing variables
  mutate(Depth = (Lower_Z-Upper_Z)/2,
         Ship = case_when(Year < 1969 ~ "Ship1",
                          Year >= 1969 ~ "Ship2"),
         Gear = case_when(Year < 1966 ~ 101,
                         Year== 1966 | Year== 1967 ~ 268,
                         Year > 1967 ~ 269),
         ShipCruise = str_c("OSP",as.character(Year),str_pad(ShipCruise, 2, pad = "0"),sep = "_"),
         UniqueSampleID = str_c(ShipCruise,str_pad(SampleID, 4, pad = "0"),as.character(Month), as.character(Day),sep = "_"))

Papa_sum <- Papa %>% 
  group_by(UniqueSampleID, Group) %>% 
  dplyr::summarise(TotAbundance = sum(as.numeric(Abundance))) # Calculate total abundance for each UniqueSampleID

Papa <- Papa %>% 
  select(-Abundance)

Papa <- left_join(Papa, Papa_sum, by = c("UniqueSampleID", "Group")) %>% 
  distinct()

Papa_wide <- Papa %>% 
  mutate(Group = as.factor(Group),
         UniqueSampleID = as.factor(UniqueSampleID)) %>% 
  group_by(UniqueSampleID) %>% 
  spread(Group, TotAbundance, fill = 0) %>% 
  ungroup()

Papa_long <- Papa_wide %>% 
  gather(c("Chaetognaths", "Euphausiids", "Larvaceans", 
                  "Salps", "Jellyfish", "OmniCopepods", 
                  "CarnCopepods"), key= "Group", value = "TotAbundance")

# Not everything was looked for all the time, 
# for five major taxa (copepods, chaetognaths, euphausiids, amphipods, and Aglantha) from 
# 1964 to 1967, and species identifications, counts and lengths for many 
# samples collected between 1968 to 1980. 

Papa <- Papa_long %>% 
  filter((Year < 1968 & Group != "OmniCopepods" & TotAbundance > 0) | Year >= 1968) %>% 
  filter((Year < 1968 & Group != "CarniCopepods" & TotAbundance > 0) | Year >= 1968) %>% 
  filter((Year < 1968 & Group != "Salps" & TotAbundance > 0) | Year >= 1968) %>% 
  filter((Year < 1968 & Group != "Larvaceans" & TotAbundance > 0) | Year >= 1968)




# I could not get a random number to work here for some reason. 
# So I set it to 9 for the moment.
Papa <- Papa %>%
  mutate(TimeLocal = replace_na(9))

# Papa <- Papa %>% 
#   mutate(TimeLocal = replace_na(TimeLocal, runif(222, min = 6, max = 12)))
# 
# Papa <- Papa %>% 
#   mutate(TimeLocal = case_when(is.na(TimeLocal) ~ runif(sum(is.na(TimeLocal)), min = 6, max = 12),
#                                !is.na(TimeLocal) ~ TimeLocal))
       

Papa <- Papa %>% 
  mutate(Group = as.factor(Group),
         ShipCruise = as.factor(ShipCruise),
         Transect = as.factor(Transect),
         UniqueSampleID = as.factor(UniqueSampleID),
         Project = as.factor(Project),
         DBase = as.factor(DBase),
         Group = as.factor(Group)) %>%
  dplyr::select(Retain) %>% 
  droplevels()

saveRDS(Papa,'OSPapa_Final.rds')



           L4 <- L4 %>% 
           
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
         saveRDS(L4long,'DatabaseOutput/L4_Final.rds')
         
         
         