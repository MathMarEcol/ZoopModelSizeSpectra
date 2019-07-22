fZerosByUniqueCruise <- function(dat,uni){
  
  for (i in 1:length(uni)){
    
    dat2 <- dat %>% 
      filter(ShipCruise == uni[i]) %>%
      ungroup()
    
    wide <- dat2 %>% spread(Group,TotAbundance, fill = 0) %>% # Wide format: Move the Abundance into each group. 
      arrange(ShipCruise, Year, Month, Day, TimeGMT) # Arrange to make it easier to debug
    
    # Calculate total abundance for each UniqueSampleID
    dfcols <- names(wide)
    wide2 <- wide %>% 
      group_by(UniqueSampleID) %>% # Regroup by UniqueSampleID
      summarize(Chaetognaths = ifelse("Chaetognaths" %in% names(wide), sum(Chaetognaths), NA),
                CarnCopepods = ifelse("CarnCopepods" %in% names(wide), sum(CarnCopepods), NA),
                OmniCopepods = ifelse("OmniCopepods" %in% names(wide), sum(OmniCopepods), NA),
                Euphausiids = ifelse("Euphausiids" %in% names(wide), sum(Euphausiids), NA),
                Jellyfish = ifelse("Jellyfish" %in% names(wide), sum(Jellyfish), NA),
                Larvaceans = ifelse("Larvaceans" %in% names(wide), sum(Larvaceans), NA),
                Salps = ifelse("Salps" %in% names(wide), sum(Salps), NA))
    
    # Remove the Groups from CPR_wide so we don't get duplicates on the inner_join
    if("Chaetognaths" %in% names(wide)) wide <- dplyr::select(wide,-c("Chaetognaths"))
    if("CarnCopepods" %in% names(wide)) wide <- dplyr::select(wide,-c("CarnCopepods"))
    if("OmniCopepods" %in% names(wide)) wide <- dplyr::select(wide,-c("OmniCopepods"))
    if("Euphausiids" %in% names(wide)) wide <- dplyr::select(wide,-c("Euphausiids"))
    if("Jellyfish" %in% names(wide)) wide <- dplyr::select(wide,-c("Jellyfish"))
    if("Larvaceans" %in% names(wide)) wide <- dplyr::select(wide,-c("Larvaceans"))
    if("Salps" %in% names(wide)) wide <- dplyr::select(wide,-c("Salps"))
    
    wide3 <- inner_join(wide, wide2, by="UniqueSampleID") # Join wide (metadata) with wide2 (Abbundance)
    
    long <- gather(wide3, Group, TotAbundance, Chaetognaths:Salps, factor_key = TRUE, na.rm = T) %>% 
      arrange(ShipCruise, Year, Month, Day, TimeGMT) %>% # Arrange to make it easier to debug 
      distinct(Group, UniqueSampleID, .keep_all = TRUE)
    
    if (i == 1) dat0 <- long else dat0 <- bind_rows(dat0, long)
    rm(long, wide, wide2, wide3)    
  }
  
  return(dat0)
}