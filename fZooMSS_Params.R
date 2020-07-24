## Set up Model Parameter List, this imports "Groups", but also sets parameters that are
## fixed across all Groups, or required to run the model

## Last Updated Tuesday 17th March 2020
##

fZooMSS_Params <- function(Groups, input_params){

  param <- list(
    Groups = Groups, # Read in functional group specific parameters from file
    ngrps = dim(Groups)[1], # no. of Groups
    dx = 0.1, # log10 weight step
    day = 12, # day length (hours of each day in sun)
    tmax = input_params$tmax, # max years
    dt = input_params$dt, # timestep
    w0 = 10^(min(Groups$W0)),		# minimum dynamic size class
    wMax = 10^(max(Groups$Wmax)),# maximum dynamic size class
    gge_base = 0.25, # baseline gross growth efficiency
    ZSpre = 1, # senescence mortality prefactor
    ZSexp = 0.3, # senescence mortality exponent
    f_mort = 0, # fishing mortality (yr^-1)
    w0_phyto = 10^(-14.5), # minimum phytoplankton size class (1um)
    wMax_phyto = 10^input_params$phyto_max, # maximum phytoplankton size class
    zoo_grps = which(Groups$Type == "Zooplankton"), # Which rows are zooplankton
    fish_grps = which(Groups$Type == "Fish"), # Which rows are fish
    num_zoo = sum(Groups$Type == "Zooplankton"), # How many zooplankton
    num_fish = sum(Groups$Type == "Fish"), # How many fish
    cc_phyto = 0.1, # Carbon content of phytoplankton size classes
    isave = 100 # how often to save results every 'isave' time steps
  )

  ## Add additional parameters which are based on the parameter set
  param2 <- list(
    w = 10^(seq(from = log10(param$w0), to = log10(param$wMax), param$dx)), # Set up dynamic weight grid
    w_phyto = 10^(seq(from = log10(param$w0_phyto), to = log10(param$wMax_phyto), param$dx)), # Set up phytoplankton size classes
    nsave  = floor(param$tmax/(param$dt*param$isave)) # Number of time slots to save
  )

  param2$ngrid <- length(param2$w) # total number of size classes for zoo and fish
  param2$ngridPP <- length(param2$w_phyto) # total number of size classes for phyto

  param <- c(input_params, param, param2) # Join with input_params
  return(param)
}
