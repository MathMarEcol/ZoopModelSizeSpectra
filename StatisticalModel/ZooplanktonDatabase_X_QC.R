
COPEPOD <- readRDS("COPEPOD_Final.rds")

dat <- COPEPOD %>% filter(Tow == "V") %>% 
  filter(Group=="Euphausiids") %>% 
  droplevels()
  
  dat2 <- dat %>% 
    group_by(Project) %>% 
    summarise(n = n(),
              zer = sum(TotAbundance==0),
              prop = (sum(TotAbundance==0))/n()) %>% 
    arrange(prop)
    
  
  dat2 <- COPEPOD %>% 
    filter(Group=="Euphausiids") %>% 
    filter(Project == "GILL Zooplankton Collection")

unique(dat$Gear)
dat$Gear <- as.factor(dat$Gear)
table(dat$Gear)

ggplot(data=dat,aes(x=1:length(TotAbundance),y=TotAbundance)) +
  geom_point()

             

