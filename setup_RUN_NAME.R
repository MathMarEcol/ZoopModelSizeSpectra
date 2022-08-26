## This is the development version of the model published in Heneghan et al., (2020):
## This version models multiple zooplankton functional groups, and three fish groups
##
## If you save timesteps, you can examine the model dynamics at
##   https://jaseeverett.shinyapps.io/ZooMSS_Dashboard/
##  Alternatively you can clone the ZooMSS_Dashboard from Github
##   https://github.com/MathMarEcol/ZooMSS_Dashboard
##
##
## Code written by Dr Jason Everett (UQ/UNSW/CSIRO), Dr Ryan Heneghan (QUT) and Mr Patrick Sykes (UQ)
## Last updated 7th October 2021


# library(Rcpp) # Only needed if we are running with Rcpp code.
source("fZooMSS_Model.R") #source the model code
source("fZooMSS_Xtras.R") #source the helper functions

# enviro_data <- readRDS("envirodata_fiveDeg_20200317.rds") # Load environmental data.

# You can also create your own environmental data using the below.
enviro_data <- fZooMSS_CalculatePhytoParam(data.frame(cellID = c(1, 2, 3), # Increment cellID so the savename changes
                                                      sst = c(5, 5, 5), # sea surface temperature
                                                      chlo = c(0.1, 0.5, 2), # Chlorophyll-a concentration
                                                      dt = c(0.01, 0.01, 0.01))) # time step to use (in years)

enviro_data$tmax <- 100 # Set length of simulation (years)

jobname <- "DATE_JOBNAME"  # This is the job name used on the HPC queue, and also to save the run: Recommend: YYYYMMDD_AbbrevExperimentName.
HPC <- FALSE # Is this being run on a HPC for all cells or will we manually choose the row of the enviro_data to be used.
SaveTimeSteps <- TRUE # Should we save all time steps. This can be very large if tmax is large
Groups <- read.csv("TestGroups_Feeding.csv", stringsAsFactors = FALSE) # Load in functional group information. This can be edited directly.
# Groups <- read.csv("TestGroups.csv", stringsAsFactors = FALSE) # Load in functional group information. This can be edited directly.


enviro_row <- 1 # Which row of the environmental data do you want to run if HPC=FALSE.

### No need to change anything below here.
if (HPC == TRUE){
  ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX')) # Get the array run number on HPC
} else {
  ID <- enviro_row
}
ID_char <- sprintf("%04d",ID) # Set the ID as a 4 digit character so it will sort properly

input_params <- enviro_data[ID,] # pick out the environmental data to use in the model

out$model$model_runtime <- system.time(  # saves the time taken to run the model
  out <- fZooMSS_Model(input_params, Groups, SaveTimeSteps) # run the model with selected parameters
)

# Save the output if you want
saveRDS(out, file = paste0("RawOutput/", jobname, "_", ID_char,".RDS"))


## Plotting
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(patchwork)
source("fZooMSS_Plot.R")

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
