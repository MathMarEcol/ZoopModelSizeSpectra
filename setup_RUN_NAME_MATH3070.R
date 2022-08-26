## This is the development version of the model published in Heneghan et al., (2020):
## This version models multiple zooplankton functional groups, and three fish groups
## Slightly modified for MATH3070
##
## If you save timesteps, you can examine the model dynamics at
##   https://jaseeverett.shinyapps.io/ZooMSS_Dashboard/
##  Alternatively you can clone the ZooMSS_Dashboard from Github
##   https://github.com/MathMarEcol/ZooMSS_Dashboard
##
##
## Code written by Dr Jason Everett (UQ/UNSW/CSIRO), Dr Ryan Heneghan (QUT) and Mr Patrick Sykes (UQ)
## Last updated 24 August 2022

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(patchwork)

source("fZooMSS_Plot.R") # Source plotting functions
source("fZooMSS_Model.R") # Source the model code
source("fZooMSS_Xtras.R") # Source the helper functions

# Set up environmental data
enviro_data <- fZooMSS_CalculatePhytoParam(data.frame(cellID = 1, # Increment cellID so the savename changes
                                                      sst = 15, # sea surface temperature
                                                      chlo = 0.2, # Chlorophyll-a concentration
                                                      dt = 0.01)) # time step to use (in years)

SaveTimeSteps <- TRUE # Should we save all time steps. This can be very large if tmax is large
enviro_data$tmax <- 100 # Set length of simulation (years)

# Read in the parameter values
Groups <- read.csv("TestGroups.csv", stringsAsFactors = FALSE) # Load in functional group information. This can be edited directly.

jobname <- "Default"  # This is the job name used on the HPC queue, and also to save the run: Recommend: YYYYMMDD_AbbrevExperimentName.

out$model$model_runtime <- system.time(                     # Saves the time taken to run the model
  out <- fZooMSS_Model(enviro_data, Groups, SaveTimeSteps)  # Run the model with selected parameters
)

# Save the output if you want
saveRDS(out, file = paste0("RawOutput/", jobname, ".RDS"))


## Plotting

# Plot predator-prey
# Have a look at the predator-prey interactions
(ggPP <- fZooMSS_PlotPredPrey(Groups))

# Plot Size Spectra
#(ggSizeSpec <- fZooMSS_Plot_SizeSpectra(out)) # abundance size spectra
(ggBiomassSizeSpec <- fZooMSS_Plot_BiomassSizeSpectra(out)) # biomass size spectra

# Plot PPMRs
(ggPPMR <- fZooMSS_Plot_PPMR(out))

## If you have saved the time steps you can plot the time series
ggBiomassTS <- fZooMSS_Plot_BiomassTimeSeries(out)  # biomass time series
ggGrowthTS <- fZooMSS_Plot_GrowthTimeSeries(out) # growth rate time series
#ggPredTS <- fZooMSS_Plot_PredTimeSeries(out) # predation mortality time series
(ggBiomassTS / ggGrowthTS) + plot_layout(guides = "collect") # plot them on top of each other

# Get abundance by size class
Ssize <- fZooMSS_SumSpecies(list(out$abundances))[[1]] # Function returns a list so we get the first

# Get abundance by species
Sspecies <- fZooMSS_SumSize(list(out$abundances))[[1]] # Function returns a list so we get the first

# Get the wet weight of each size class
WWsize <- fZooMSS_SizeBiomass(list(out$abundances), out$model)[[1]] # Function returns a list so we get the first

# Get the wet weight of each group
WWspecies <- fZooMSS_SpeciesBiomass(list(out$abundances), out$model)[[1]] # Function returns a list so we get the first

# Get the carbon biomass of each group
C <- fZooMSS_SpeciesCarbonBiomass(list(out$abundances), out$model)[[1]] # Function returns a list so we get the first

# Get the trophic level of each group
TL <- fZooMSS_TrophicLevel(out)

df <- data.frame(Groups, Abundance = Sspecies, Biomass = WWspecies, TrophicLevel = TL, CarbonBiomass = C)
write.csv(x = df, file = paste0("Outputs", jobname, ".csv"))

sum(WWspecies[3:9])
sum(WWspecies[10:12])
