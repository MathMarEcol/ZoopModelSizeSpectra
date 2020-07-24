# This script contains a list of helper functions which can be used to analyse and plot the ZooMSS output

# Summarise the biomass for each gridcell by species
fZooMSS_SpeciesBiomass = function(res,w) {
  Biomass <- map(res, function(x) apply(sweep(x, 2, w, '*'), 1, sum))
  return(Biomass)
}

# Convert Abundance to Biomass for all species and weight classes
fZooMSS_Biomass <- function(res, w){
  Biomass <- map(res, function(x) sweep(x, 2, w, '*')) # Biomass in grams
  return(Biomass)
}

# Summarise the biomass for each gridcell by sizeclass
fZooMSS_SizeBiomass = function(res,w) {
  Biomass <- map(res, function(x) apply(sweep(x, 2, w, '*'), 2, sum))
  return(Biomass)
}



# Function to calculate the mean of the last 50 % of the model
fZooMSS_AveOutput = function(x){
  ave_x <- colMeans(x[(ceiling(0.5*(dim(x)[1])):dim(x)[1]),,], dims = 1)
  return(ave_x)
}

# Remove nonsense attributes if we are working for speed and memory efficiency.
untibble <- function (tibble) {
  data.frame(unclass(tibble), check.names = FALSE, stringsAsFactors = FALSE)
}  ## escape the nonsense





