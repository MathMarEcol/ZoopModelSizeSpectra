### Run the multi-zoo model
rm(list = ls())

setwd("~/Desktop/Multi-Zoo Model") ## CHANGE YOUR WORKING DIRECTORY
source("Multi-Zoo Model.R")

enviro_data <- readRDS("store5.RDS")
enviro_data$dt <- 0.01

#a <- Sys.time()
Groups <- read.csv("Test Groups.csv")
enviro <- enviro_data[441,]
param <- params(Groups, enviro, tmax = 300, f_mort = 0) # Set up parameter list
model <- Setup(param) # Set up model equation stuff 
modelss <- Project(model, fish_on = TRUE) # Run the model
#b <- Sys.time()
#b-a

Spectrum_Plot(modelss, param, fish_on = TRUE)
Biomass_Plot(modelss, param)

Bio_Cont_Plot(modelss, param)
#Fish_Diet(modelss, param)
PPMR_plot(modelss, param)

general_results <- Summary_Results(modelss, param)$General
zoo_results <- Summary_Results(modelss, param)$Zoo_Groups

### MAKE A PIECHART OF COMPOSITIONS
# -4.8 is 300um, so biomass/abundance above 200um
Zoo_Biom_Pie(modelss, -8, 0) 
Zoo_Abund_Pie(modelss, -8, 0)
