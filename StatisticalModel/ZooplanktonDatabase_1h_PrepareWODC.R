rm(list = ls())

library(tidyverse)
library(lubridate)

Retain <- c("ShipCruise","Transect","UniqueSampleID","Year","Month","Day","TimeGMT","TimeLocal",
            "Latitude","Longitude","TotAbundance","Depth","Upper_Z","Lower_Z","Tow","Gear",
            "Type","Mesh","Volume","Project","Ship","Group","DBase")

dat <- read_delim("Data/WODC/WODC_ocldb1559734053.7056.OSD.csv", delim = "|",col_names = "Temp1")

x <- which(str_detect(dat$Temp1, "#----------------")) # Find the start of all tows

WODC <- data.frame("ShipCruise" = NA)

for (i in 1:length(x)){
  
  x1 <- x[1]+1 # Start Row
  x2 <- x[i+1]-2 # End Row
  temp <- dat[x1:x2,] # Exract each tow
  write_delim(temp,"temp.csv", delim = "|", col_names = FALSE) # Save the tow
  
  # Reload the tow as csv
  temp <- read_csv("temp.csv",col_names = FALSE)
  
  WODC$ShipCruise[i] <- temp$X3[2]
  WODC$Latitude[i] <- temp$X3[5]
  WODC$Longitude[i] <- temp$X3[6]
  
  WODC$Year[i] <- temp$X3[7]
  WODC$Month[i] <- temp$X3[8]
  WODC$Day[i] <- temp$X3[9]
  
  
  
  
}






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
saveRDS(dat,'CalCOFI_Final.rds')





