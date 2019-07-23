## An extension of the model published in Heneghan et al., (2016):
## Models multiple zooplankton functional groups, and three fish groups
## This code is to run the model across multiple cores

rm(list = ls())
library(rslurm)
setwd("~/Desktop/Multi Core Code")

enviro_data <- read.csv("enviro_5d_3.csv") # Load in environmental data
tt = which(enviro_data$sst < 1) # remove squares where sst < 1 (poor coverage)
enviro_data = enviro_data[-tt,]

lat_lon = enviro_data[,c("x","y")] # Extract lats and longs
enviro_data$dt <- 0.01 # Set time step (years) for model runs
enviro_data <- enviro_data[,c( "sst", "chlo","a","b", "phyto_max", "dt" )]
enviro_data$tmaxx <- 500 # Set length of simulation (years)

### WE WANT 500 OBSERVATIONS, WITH ALL OBSERVATIONS ABOVE 1MG M-3 CHLO
#e_large = enviro_data[which(log10(enviro_data$chlo) >= 0),]
#e_all_else = enviro_data[-which(log10(enviro_data$chlo) >= 0),]
#e_all_else_samp = e_all_else[sample(seq(1,dim(e_all_else)[1]),400, replace = FALSE),]
#enviro_data_samp = rbind(e_all_else_samp, e_large)
#saveRDS(enviro_data_samp, file = "store5.RDS")

#enviro_data <- readRDS("store5.RDS") # Run this if you only want 503 observations, instead of full ~1600
#enviro_data <- enviro_data[c(1:150, 451:500),]
#enviro_data$tmaxx <- 500 # Set length of simulation (years)

Groups <- read.csv("Test Groups.csv") # Load in functional group information

## To run in slurm, the model needs to be wrapped in a function, this is "multi_zoo_slurm" here
multi_zoo_slurm <- function(sst, chlo, a, b, phyto_max, dt, tmaxx){
  
  ## Set up Model Parameter List, this imports "Groups", but also sets parameters that are
  ## fixed across all groups, or required to run the model
  params <- function(fileGroups, enviroo, tmax, f_mort){
    
    groups = fileGroups  # Read in functional group specific parameters from file
    nutrition = groups$carbon # Extract carbon content (nutr. quality) of each group
    environ = enviroo # Environmental information
    
    # Set up parameter list
    param = list(groups = groups, 
                 environ = environ,
                 nutrition = nutrition,
                 
                 # Model parameters
                 ngrps = dim(groups)[1],		# no. of groups
                 tmax = tmax,					# no. of years
                 dt = environ$dt, 				# time step
                 dx = 0.1,         # log10 weight step
                 day = 12,          # day length (hours of each day in sun)
                 gge_base = 0.25, # baseline gross growth efficiency
                 w0 = 10^(min(groups$W0)),		# minimum dynamic size class
                 wMax = 10^(max(groups$Wmax)),# maximum dynamic size class
                 ZSpre = 1, # senescence mortality prefactor
                 ZSexp = 0.3, # senescence mortality exponent
                 f_mort = f_mort, # fishing mortality (yr^-1)
                 w0_phyto = 10^(-14.5),		# minimum phytoplankton size class (1um)
                 wMax_phyto = 10^environ$phyto_max		# maximum phytoplankton size class
    )
    
    param$isave = 100	# how often to save results every 'isave' time steps
    param$fish_grps = which(is.na(groups$prop) == TRUE) # Which rows are fish
    param$zoo_grps = which(is.na(groups$prop) == FALSE) # Which rows are zooplankton
    param$num_zoo = sum(is.na(param$groups$prop) == FALSE) # How many zooplankton
    param$num_fish = sum(is.na(param$groups$prop) == TRUE) # How many fish
    return(param)
  }
  
  ## The Setup function calculates the feeding kernels, dvm (if included), temp effects...
  ## everything that can be solved before we start iterating through time to solve MvF-D
  Setup = function(param){
    
    # Pull out some useful parameters - just a shortcut
    grp = param$groups # Functional group info table
    ngrps = param$ngrps # number of functional groups
    dt = param$dt # time step
    dx = param$dx # log10 weight step
    environ = param$environ # environmental info
    fish_grps = param$fish_grps # which rows of grp are fish
    zoo_grps = param$zoo_grps # which rows of grp are zoo
    num_zoo = param$num_zoo # how many zooplankton groups
    num_fish = param$num_fish # how many fish groups
    
    # Set up dynamic weight grid
    w = 10^(seq(from = log10(param$w0), to =  log10(param$wMax), dx))
    ngrid = length(w) # total number of size classes for zoo and fish
    
    # Set up phytoplankton size classes
    w_phyto = 10^(seq(from = log10(param$w0_phyto), to = log10(param$wMax_phyto), dx))
    ngridPP = length(w_phyto) # total number of size classes for phyto
    
    # Number of time slots to save
    nsave   = floor(param$tmax/(param$dt*param$isave))
    
    
    # Diel Vertical Migration - change availability of phyto to zoo
    # and zoo to fish based on slope of phytoplankton (calculated as
    # proportion of day searching for food and available for predation)
    dvm_max = param$day/24 # maximum proportion of day spent migrating
    ESD_sizes = 2*(3/(4*pi)*w)^(1/3) # convert g wet weight to ESD (cm)
    dvm = dvm_max*(1.02*(ESD_sizes) - 0.02) # size-dependent amount of time spent away from surface
    dvm[which(w < 10^-5.4)] = 0 # Microzoo don't migrate (ESD < 0.02cm)
    dvm[which(w > 10^-0.3)] = dvm_max # Macrozoo migrate max time (ESD > 10cm)
    
    #plot(log10(w), dvm, ylab = "Proportion of time spent DVMing (away)") # If you want to have a look at how dvm is parameterised across body size
    
    dvm_mat = matrix(dvm, nrow = ngrid, ncol = ngrid, byrow = TRUE) # Matrix of dvm, nrow = number of pred size classes
                                                                    # ncol = number of pre size classes
    
    # This works out the proportion of time an predator of size w will have access to a prey of size w', for all w and w'
    dvm_mat = 1 -  dvm_mat 
    dvm_mat[lower.tri(dvm_mat)] = 0
    dvm_mat = t(dvm_mat) + dvm_mat
    diag(dvm_mat) = diag(dvm_mat)/2

    # Dynamic prey availability matrix: dim1 is predators, dim2 is predator size classes,
    # dim3 is prey groups, dim 4 is prey size classes.
    dynam_theta = array(1, dim = c(ngrps, ngrid, ngrps, ngrid))
    dynam_theta = sweep(dynam_theta, c(2,4), dvm_mat,"*")
    
    # Phyto availability matrix: rows are predators, columns are their size classes,
    # entries are time spent feeding on phytoplankton for the size class
    phyto_theta = matrix(1-dvm, nrow = ngrps, ncol = ngrid, byrow = TRUE)
    
    ### REMOVE DVM
   dynam_theta = array(1, dim = c(ngrps, ngrid, ngrps, ngrid))
   phyto_theta = matrix(1, nrow = ngrps, ncol = ngrid, byrow = TRUE)
    ###
    
    carn_grps = which(grp$type == 'C')
    phyto_theta[carn_grps,] = 0 # Carnivorous groups can't eat phyto
    
    cc_phyto = 0.1   # Carbon content of phytoplankton size classes
    
    ## Makes the model object, full of constant functions for model
    model = list(
      param = param,
      environ = environ,
      ngrid = ngrid,
      ngridPP = ngridPP,
      
      # Phytoplankton abundance
      nPP = 10^(environ$a)*(w_phyto^(environ$b)), # phyto abundance spectrum
      
      # Grid parameters
      w = w,
      dx = dx,
      w_phyto = w_phyto,
      
      # Group parameters storage
      phyto_growthkernel = array(NA, dim = c(ngrps, ngrid, ngridPP)), # predation on phytoplankton
      phyto_diffkernel = array(NA, dim = c(ngrps, ngrid, ngridPP)), # diffusion from phytoplankton consumption
      phyto_dietkernel = array(NA, dim = c(ngrps, ngrid, ngridPP)), # diet from phytoplankton
      dynam_growthkernel =  array(NA, dim = c(ngrps, ngrid, ngrid)), # predation on zoo and fish
      dynam_diffkernel = array(NA, dim = c(ngrps, ngrid, ngrid)), # diffusion from zoo and fish consumption
      dynam_dietkernel = array(NA, dim = c(ngrps, ngrid, ngrid)), # diet from zoo and fish
      dynam_mortkernel = array(NA, dim = c(ngrps, ngrid, ngrid)), # mortality from predation on dynamic component
      M_sb = matrix(0, nrow = ngrps, ncol = ngrid), # senescence mortality
      fish_mort = matrix(0, nrow = ngrps, ncol = ngrid), # fishing mortality
      
      # Output storage
      N = array(0, dim = c(nsave, ngrps, ngrid)), # dynamic abundance spectrum
      Z = array(0, dim = c(nsave, ngrps, ngrid)), # total mortality
      gg = array(0, dim = c(nsave, ngrps, ngrid)), # growth
      diet = array(0, dim = c(nsave, c(ngrps), c(ngrps+3))), # diet
      Biomass = matrix(0, nrow = nsave, ncol = ngrps), # biomass of each group
      Diff = array(0, dim = c(nsave, ngrps, ngrid)) # save diffusion
    )
    
    # GGE for different groups
    assim_phyto =  (param$groups$alpha)*cc_phyto # Phytoplankton
    assim_dynam =  matrix(param$groups$alpha*param$nutrition, nrow = ngrps, ncol = ngrps, byrow = TRUE) # rows are predators, columns are prey
    
    #### INITIAL DYNAMIC POPULATION ABUNDANCES
    a_dynam = 10^(environ$a)*(w[1]^(environ$b+1)) # calculate coefficient for initial dynamic spectrum, so that N(w_phyto) equals
                                                  # N(w_dynam) at w[1]
    
    # Initial abundances form a continuation of the plankton spectrum, with a slope of -1
    tempN = matrix(a_dynam*(w)^-1, nrow = ngrps, ncol = ngrid, byrow = TRUE) 
    props_z = grp$prop[zoo_grps] # Zooplankton proportions
    tempN[zoo_grps,] = props_z*tempN[zoo_grps,] # Set abundances of diff zoo groups based on smallest size class proportions
    tempN[fish_grps,] = (1/num_fish)*tempN[fish_grps,] # Set abundandances of fish groups based on smallest size class proportions
    
    # For each group, set densities at w > Winf and w < Wmin to 0
    tempN[unlist(tapply(round(log10(w), digits = 2), 1:length(w), function(wx,Winf) Winf < wx, Winf = (grp$Wmax)))] = 0
    tempN[unlist(tapply(round(log10(w), digits = 2), 1:length(w), function(wx,Wmin) Wmin > wx, Wmin = (grp$W0)))] = 0
    model$N[1,,] = tempN
    
    # Fishing mortality
    model$fish_mort[fish_grps, c(w >= 1)] = param$f_mort
    
    ### MATRICES FOR LOG TRANSFORM OF EQUATION
    # Predators are rows, phyto prey weights are columns   
    gg_log_t_phyto = ((w^-1) %*% t(w_phyto))/log(10) # Growth
    diff_log_t_phyto = ((w^-2) %*% t(w_phyto^2))/log(10) # Diffusion
    diet_log_t_phyto = matrix(w_phyto, nrow = length(w), ncol = length(w_phyto), byrow = TRUE) # Diet/Ingestion
    
    # Predators are rows, dynam prey weights are columns
    gg_log_t_dynam = ((w^-1) %*% t(w))/log(10) # Growth
    diff_log_t_dynam = ((w^-2) %*% t(w^2))/log(10) # Diffusion
    diet_log_t_dynam = matrix(w, nrow = length(w), ncol = length(w), byrow = TRUE) # Diet/ingestion
    
    ### PREDATION KERNELS FOR PHYTOPLANKTON SPECTRUM AND DYNAMIC SPECTRUM
    phyto_pred_weight_matrix = matrix(w, nrow = ngrid, ncol = ngridPP)
    dynam_pred_weight_matrix = matrix(w, nrow = ngrid, ncol = ngrid)
    phyto_prey_weight_matrix = matrix(w_phyto, nrow = ngrid, ncol = ngridPP, byrow = TRUE)
    dynam_prey_weight_matrix = matrix(w, nrow = ngrid, ncol = ngrid, byrow = TRUE)
    
    ## Search Volume storage
    SearchVol = matrix(NA, nrow = ngrps, ncol = ngrid) # Search volume
    
    # Simpson's Rule matrices for growth, diffusion and mortality integrals
    simp_phyto = array(1, dim = ngridPP)
    simp_phyto[c(seq(2,ngridPP-1,2))] = 4
    simp_phyto[c(seq(3,ngridPP-1,2))] = 2
    sm_phyto = matrix(simp_phyto, nrow = ngrid, ncol = ngridPP, byrow = TRUE)*(dx/3)
    
    simp_dynam = array(1, dim = ngrid)
    simp_dynam[c(seq(2,ngrid-1,2))] = 4
    simp_dynam[c(seq(3,ngrid-1,2))] = 2
    sm_dynam = matrix(simp_dynam, nrow = ngrid, ncol = ngrid, byrow = TRUE)*(dx/3)
    
    ## Temperature Effect Matrix
    # Effect of temperature on feeding and predation rate
    #temp_flag <- 2.4^((environ$sst)/10)
    #temp_cil <-  2.8^((environ$sst)/10)
    #temp_gel <- 1.75^((environ$sst)/10)
    #temp_larv <- 2.2^((environ$sst)/10)
    #temp_chaet <- 2.44^((environ$sst)/10)
    #temp_crust <- 2.57^((environ$sst)/10)
    #temp_ocop <- 2.8^((environ$sst)/10)
    #temp_ccop <- 2.8^((environ$sst)/10)
    #temp_fish <- 2.6^((environ$sst)/10)
    
    #temp_zoo <- c(temp_flag, temp_cil, temp_larv, temp_ocop, temp_ccop, temp_crust, temp_chaet, temp_larv, temp_gel)
    #temp_fish <- rep(temp_fish, num_fish)
    
    #temp_zoo <- rep(exp(23.93 - 0.59/(8.62e-05*(273+environ$sst))), num_zoo) # exp(23.93 - 0.59/(8.62e-05*(273+environ$sst)))
    #temp_fish <- rep(exp(23.93 - 0.59/(8.62e-05*(273+environ$sst))), num_fish) # exp(25.55 - 0.63/(8.62e-05*(273+environ$sst)))/exp(23.93 - 0.59/(8.62e-05*(273+environ$sst)))
    
    ### Q10 OF 2 FOR ALL ZOO AND FISH
    temp_zoo <- rep(2.^((environ$sst - 30)/10), num_zoo) # exp(23.93 - 0.59/(8.62e-05*(273+environ$sst)))
    temp_fish <- rep(2.^((environ$sst - 30)/10), num_fish)
    temp_effect <- matrix(c(temp_zoo, temp_fish), nrow = ngrps, ncol = ngrid)
    
    #### CALCULATES CONSTANT BITS OF THE MODEL FUNCTIONS FOR EACH GROUP
    for(i in 1:ngrps){ 
      ## Senescence mortality
      if(i < 10){
      model$M_sb[i,] = param$ZSpre*(w/(10^(grp$Wmat[i])))^param$ZSexp
      model$M_sb[i, 10^(grp$Wmax[i]) < w] = 0
      model$M_sb[i, 10^(grp$Wmat[i]) > w] = 0	
      }
      
      if(i > 9){
      model$M_sb[i,] = 0.1*param$ZSpre*(w/(10^(grp$Wmat[i])))^param$ZSexp
      model$M_sb[i, 10^(grp$Wmax[i]) < w] = 0
      model$M_sb[i, 10^(grp$Wmat[i]) > w] = 0	  
      }
      
      ### Search volume
      SearchVol[i,] = (grp$gamma[i])*(w^(grp$q[i]))
      SearchVol[i, 10^(grp$Wmax[i]) < w] = 0
      SearchVol[i, 10^(grp$W0[i]) > w] = 0		
      
      ### Predation Kernels
      if(is.na(grp$m[i]) == FALSE){ # If group has an m-value (zooplankton)
        # Calculate PPMR for zooplankton, which changes according to body-size (Wirtz, 2012)	
        D.z = 2*(3*w*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
        betas =  (exp(0.02*log(D.z)^2 - grp$m[i] + 1.832))^3 # Wirtz's equation
        beta_mat_phyto = matrix(betas, nrow = ngrid, ncol = ngridPP) 
        beta_mat_dynam = matrix(betas, nrow = ngrid, ncol = ngrid)
        
        # Calculate feeding kernels
        sp_phyto_predkernel = exp(-0.5*(log((beta_mat_phyto*phyto_prey_weight_matrix)/
                                              phyto_pred_weight_matrix)/grp$sigma[i])^2)/
          sqrt(2*pi*grp$sigma[i]^2)
        sp_dynam_predkernel = exp(-0.5*(log((beta_mat_dynam*dynam_prey_weight_matrix)/
                                              dynam_pred_weight_matrix)/grp$sigma[i])^2)/
          sqrt(2*pi*grp$sigma[i]^2)
        
      } else { # If group does not have an m-value (fish)
        beta_mat_phyto = matrix(grp$beta[i], nrow = ngrid, ncol = ngridPP)
        beta_mat_dynam = matrix(grp$beta[i], nrow = ngrid, ncol = ngrid)
        
        # Calculate feeding kernels
        sp_phyto_predkernel = exp(-0.5*(log((beta_mat_phyto*phyto_prey_weight_matrix)/
                                              phyto_pred_weight_matrix)/grp$sigma[i])^2)/
          sqrt(2*pi*grp$sigma[i]^2)
        sp_dynam_predkernel = exp(-0.5*(log((beta_mat_dynam*dynam_prey_weight_matrix)/
                                              dynam_pred_weight_matrix)/grp$sigma[i])^2)/
          sqrt(2*pi*grp$sigma[i]^2)
      }
      
      ### GROWTH INTEGRAL CONSTANTS
      # Predators are rows, prey are columns	
      model$phyto_growthkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP)*
        sp_phyto_predkernel*gg_log_t_phyto*sm_phyto 
      model$dynam_growthkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
        sp_dynam_predkernel*gg_log_t_dynam*sm_dynam
      
      ### DIET INTEGRAL CONSTANTS
      # Predators are rows, prey are columns
      model$phyto_dietkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP)*
        sp_phyto_predkernel*diet_log_t_phyto*sm_phyto
      model$dynam_dietkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
        sp_dynam_predkernel*diet_log_t_dynam*sm_dynam
      
      ### DIFFUSION INTEGRAL CONSTANTS
      # Predators are rows, prey are columns	
      model$phyto_diffkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngridPP)*
        sp_phyto_predkernel*diff_log_t_phyto*sm_phyto
      model$dynam_diffkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid)*
        sp_dynam_predkernel*diff_log_t_dynam*sm_dynam	
      
      ### MORTALITY INTEGRAL CONSTANTS                                  
      # Prey are rows, predators are columns
      model$dynam_mortkernel[i,,] = matrix(SearchVol[i,], nrow = ngrid, ncol = ngrid, byrow = TRUE)*
        t(sp_dynam_predkernel)*sm_dynam 
    }
    
    no_sen = which(grp$species == c("Flagellates", "Ciliates")) # no senescence mortality for flagellates and ciliates
    #model$M_sb[c(ngrps)] = 0 
    model$M_sb = temp_effect*model$M_sb # Incorporate temp effect on senscence mortality
    
    ## Incorporate dvm, temperature effects and gross growth efficiency (assim)
    model$phyto_growthkernel = sweep(sweep(model$phyto_growthkernel, c(1,2), phyto_theta, "*"), 1, assim_phyto, "*")
    model$phyto_diffkernel = sweep(sweep(model$phyto_diffkernel, c(1,2), phyto_theta, "*"), 1, assim_phyto^2, "*")
    model$phyto_dietkernel =  sweep(sweep(model$phyto_dietkernel, c(1,2), phyto_theta, "*"), 1, 1, "*")
    
    # Dim 1 = pred group, dim2 = pred sizes, dim 3 = prey group, dim 4 = prey sizes
    model$dynam_growthkernel = sweep(sweep(sweep(dynam_theta, c(1,2,4), model$dynam_growthkernel, "*"), 
                                           c(1,3), assim_dynam, "*"), c(1,2), temp_effect, "*")
    model$dynam_diffkernel = sweep(sweep(sweep(dynam_theta, c(1,2,4), model$dynam_diffkernel, "*"), c(1,3), assim_dynam^2, "*"),
                                   c(1,2), temp_effect^2, "*")     
    model$dynam_dietkernel = sweep(sweep(sweep(dynam_theta, c(1,2,4), model$dynam_dietkernel, "*"), 
                                         c(1,3), 1, "*"), c(1,2), temp_effect, "*")
    
    # Dim 1 = prey group, dime 2 = prey sizes, dim 3 = pred group, dim 4 = pred sizes
    model$dynam_mortkernel = sweep(aperm(sweep(aperm(dynam_theta, c(3,1,4,2)), 
                                               c(2,3,4), model$dynam_mortkernel, "*"), c(1,3,2,4)),
                                   c(3,4), temp_effect, "*")
    
    #### Because phyto spectrum is constant, we can solve the phyto component of growth, and diffusion before time loop
    model$ingested_phyto = temp_effect*(rowSums(sweep(model$phyto_growthkernel, 3, model$nPP, "*"), dims = 2)) # Ingested phyto
    model$diff_phyto = temp_effect^2*(rowSums(sweep(model$phyto_diffkernel, 3, model$nPP, "*"), dims = 2)) # Diffusion from phyto
    model$diet_phyto = temp_effect*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP, "*"), dims = 2)) # Diet of total phyto
    
    ## Diet of phyto from pico, nano and micro size classes
    model$diet_pico_phyto = temp_effect*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$w_phyto) < -11.5), "*"), dims = 2))
    model$diet_nano_phyto = temp_effect*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$w_phyto) >= -11.5 & log10(model$w_phyto) < -8.5), "*"), dims = 2))
    model$diet_micro_phyto = temp_effect*(rowSums(sweep(model$phyto_dietkernel, 3, model$nPP*c(log10(model$w_phyto) >= -8.5), "*"), dims = 2))
    
    return(model)
  } # End of Setup function
  
  
  ## Run model forward in time
  Project <- function(model, fish_on){
    
    # Pull out some useful parameters - just a shortcut
    param <- model$param
    grp <- param$group
    ngrid <- model$ngrid
    ngridPP <- model$ngridPP
    ngrps <- param$ngrps
    dt <- param$dt
    fish_grps <- param$fish_grps
    zoo_grps <- param$zoo_grps
    dx <- model$dx
    w <- model$w
    w_phyto <- model$w_phyto
    dw <- model$dw
    assims <- param$nutrition/0.1*(grp$alpha)
    w0idx <- which(grp$W0 > min(grp$W0) & is.na(grp$prop) == FALSE)
    w0mins <- rep(0, length(w0idx))
    
    props_z <- grp$prop[w0idx] # Zooplankton proportions
    
    for(i in 1:length(w0idx)){ # Which size class is the smallest size class for each functional group
      w0mins[i] <- which(round(log10(w), digits = 2) == grp$W0[w0idx[i]])
    }
    
    # Handy stuff
    idx <- 2:ngrid # size class sequence
    itimemax  <- param$tmax / dt  #max index of time array
    
    # Matrices for MvF and MvF-D numeric solution
    A.iter <- C.iter <- S.iter <- A <- B <- C <- S <- matrix(0,nrow=ngrps,ncol=ngrid)
    
    # Temporary Matrices that get updated each time step
    # some of these saved for output
    N <- matrix(model$N[1,,], nrow = ngrps, ncol = ngrid) # Abundances of functional groups, dim 1 = groups, dim 2 = size classes
    nPP <- model$nPP # Abundances of phytoplankton spectrum
    
    pb = txtProgressBar(min = 0, max = itimemax, initial = 1, style = 3) # Initial progress bar
    
    # BIG TIME LOOP
    for (itime in 1:itimemax){
      
      setTxtProgressBar(pb, itime) # Update progress bar
      
      ### Create an ngrps*ngrid*ngrps*ngrid array of abundances, to save time without sweeps
      # dim1 = pred groups, dim 2 = pred sizes, dim 3 = prey groups, dim 4 = prey sizes
      N_array <- aperm(replicate(ngrid, N), c(3,1,2))
      N_array <- aperm(replicate(ngrps, N_array), c(4,1,2,3))

      
      ### GROWTH
      gg <- (model$ingested_phyto + 
               rowSums(rowSums(model$dynam_growthkernel*N_array, dims = 3), dims = 2))     
      
      ### MORTALITY
      # Predation mortality
      M2 <- (rowSums(rowSums(model$dynam_mortkernel*N_array, dims = 3), dims = 2))
      
      # Total dynamic spectrum mortality
      Z = M2 + model$M_sb  + model$fish_mort
      
      
      ### DIFFUSION
      diff <- (model$diff_phyto + rowSums(rowSums(model$dynam_diffkernel*N_array, dims = 3), dims = 2))
      
      ### MvF WITH DIFFUSION ALGORITHM
      
      idx.iter <- 2:ngrid
      idx <- 2:(ngrid-1)
      
      # Numerical implementation matrices (for MvF without diffusion)
      A.iter[,idx.iter] <- dt/dx*gg[,idx.iter-1] 
      C.iter[,idx.iter] <- 1 + dt*Z[,idx.iter] + dt/dx*gg[,idx.iter]
      S.iter[,idx.iter] <- N[,idx.iter]
      N.iter <- N
      
      # Numerical implementation matrices (for MvF WITH diffusion)
      A[,idx] <- dt/dx*(gg[,idx-1] + diff[,idx-1]*(log(10)/2+1/(2*dx)))
      B[,idx] <- diff[,idx+1]*dt/(2*dx^2)
      C[,idx] <- 1 + dt*Z[,idx] + dt/dx*(gg[,idx] + diff[,idx]*(log(10)/2+1/dx))
      S[,idx] <- N[,idx]
      
      for(i in 1:ngrps){
        ## Set size range index for current group
        curr_min_size = which(round(log10(w), digits = 2) == param$groups$W0[i])
        curr_max_size = which(round(log10(w), digits = 2) == param$groups$Wmax[i])
        idx_curr = (curr_min_size+1):curr_max_size
        
        for(j in idx_curr){## Find the abundance at the next size class with standard MvF
          N.iter[i,j] <- (S.iter[i,j] + A.iter[i,j]*N[i,j-1])/(C.iter[i,j])
          if(j >= (idx_curr[1]+1)){ ## Find abundance with MvF with diffusion
            k = j - 1
            N[i,k] = (S[i,k] + A[i,k]*N[i,k-1] + B[i,k]*N.iter[i,k+1])/C[i,k]
          }
          # MvF without diffusion for last size class  
          if(j == idx_curr[length(idx_curr)]){ 
            N[i,j] = 0
            # N[i,curr_min_size] <- N.iter[i,curr_min_size] # Keep starting sizes constant
          }
        }
      }
      
      
      #### Keep smallest fish community size class as equal to equivalent zooplankton size class
      
      ### Keep smallest zooplankton size class abundnace 
      ### for each group locked to others in size spectrum
      for(i in 1:length(w0idx)){
        w_min_curr = w0mins[i]
        exclude_mins = w0idx[which(w0mins == w_min_curr)]
        N[w0idx[i], w_min_curr] = props_z[i]*sum(N[-exclude_mins, w_min_curr])
      }
      
      
      fish_mins = unlist(lapply(param$groups$W0[fish_grps], 
                                function(x){which(round(log10(model$w), digits = 2) == x)}))
      
      if(length(fish_grps) > 1){
        N[fish_grps,fish_mins] = (1/3)*(colSums(N[-fish_grps,fish_mins]))
      }else{
        N[fish_grps, fish_mins] = sum(N[-fish_grps, fish_mins])
      }
      
      if(fish_on == FALSE){
        N[fish_grps,] = 0 # switch off fish_groups
      }
      
      # Save results:
      if((itime %% param$isave) == 0){
        isav=itime/param$isave
        
        ## Phytoplankton diet
        pico_phyto_diet = rowSums(model$diet_pico_phyto*N) # Pico-phytoplankton
        nano_phyto_diet = rowSums(model$diet_nano_phyto*N) # Nano-phytoplankton
        micro_phyto_diet = rowSums(model$diet_micro_phyto*N) # Micro-phytoplankton
        
        phyto_diet = cbind(pico_phyto_diet, nano_phyto_diet, micro_phyto_diet)
        
        ## Functional group diet
        dynam_diet =  rowSums(aperm(rowSums(sweep(model$dynam_dietkernel*N_array, c(1,2), N, "*"), dims = 3), c(1,3,2)), dims = 2)
        
        model$diet[isav,,1:3] = phyto_diet
        model$diet[isav,,c(4:15)] = dynam_diet
        
        model$N[isav,,] <- N # Save abundance
        
        ## NEED TO CHECK WHAT THIS IF STATEMENT IS FOR, IT DOESN'T APPEAR TO DO ANYTHING
        if(length(zoo_grps) > 1){
          model$N[isav,c(1,2),c(60,61)] <- 0
        }
        
        ## Save biomass
      #  model$Biomass[isav,] <- rowSums(model$N[isav,,] # Save biomass
      #                                  *matrix(model$w, nrow = ngrps, ncol = ngrid, byrow = TRUE)) 
        
        ## Save mortality rates
        #	model$Z[isav,,] <- M2 # Save total predation mortality rates
        
        ## Save growth
      	model$gg[isav,,] <-  (model$diet_phyto + 
      	                     rowSums(rowSums(model$dynam_dietkernel*N_array, dims = 3), dims = 2))     
      }
      
    } # End of time loop
    
    return(model)
    
  }
  
  
  enviro_vector <- data.frame("sst" = sst, "chlo" = chlo, "a" = a, "b" = b,
                              "phyto_max" = phyto_max, "dt" = dt)
  
  ##################### RUN THE MODEL ################################################
  param <- params(Groups, enviro_vector, tmax = tmaxx, f_mort = 0) # Set up parameter list
  model <- Setup(param) # Set up model equation stuff 
  modelss <- Project(model, fish_on = TRUE) # Run the model
  
  ################### OUTPUT ABUNDANCES ##############################################
  ave_abundances = colMeans(modelss$N[(ceiling(0.5*dim(modelss$N)[1])):(dim(modelss$N)[1]),,], dim = 1)
  ave_diets = colMeans(modelss$diet[(ceiling(0.5*dim(modelss$diet)[1])):(dim(modelss$diet)[1]),,], dim = 1)
  ave_growth = colMeans(modelss$gg[(ceiling(0.5*dim(modelss$gg)[1])):(dim(modelss$gg)[1]),,], dim = 1)
  ave_pred = colMeans(modelss$Z[(ceiling(0.5*dim(modelss$Z)[1])):(dim(modelss$Z)[1]),,], dim = 1)  
  results = list("abundances" = ave_abundances, "diets" = ave_diets, "growth" = ave_growth) #, "diets" = ave_diets
  return(results)
}

sjob <- slurm_apply(f = multi_zoo_slurm, params = enviro_data, jobname = 'multi_zoo_model',
                    add_objects = c("Groups"), nodes = dim(enviro_data)[1], cpus_per_node = 1,
                    slurm_options = list(ntasks = as.character(dim(enviro_data)[1]),
                       "mem-per-cpu" = "30000", time = "0-5:00"),
                    libPaths = "~/R_libs", submit = FALSE)

## 22/11/18 02 - Ccops from 0.4 to 0.3
## 22/11/18 02 - Larvs from 0.8 to 0.7, salps from 0.05 to 0.04, Omnicops from 0.08 to 0.07
## 23/11/18 01 - Larvs back to 0.8, salps back to 0.05, ccops to 0.25 from 0.3
## 23/11/18 02 - Time = 500 years from 1000 years, 
## 28/11/18 01 - Salps from 0.05 to 0.02, omnicops from 0.07 to 0.05
## 29/11/18 01 - Salps back to 0.05, omnicops to 0.06 from 0.07, ccops from 0.25 to 0.1
## 29/11/18 02 - Ocops from 0.07 to 0.05, larvaceans from 0.8 to 0.5
## 30/11/18 01 - Jellyfish from 0.01 to 0.02, ccops from 0.1 to 0.15, euphausiids from 0.8 to 0.7
## 30/11/18 02 - Jellyfish from 0.02 to 0.03, ccops from 0.15 to 0.13, salps from 0.05 to 0.03
## 30/11/18 03 - Krill from 0.7 to 0.2
## 30/11/18 04 - Jellyfish from 0.03 to 0.01, occops from 0.05 to 0.02
## 01/12/18 01 - Jellyfish from 0.01 to 0.03
## 11/12/18 01 - Jellyfish 0.01 from 0.03, salps 0.01 from 0.03, euphs 0.1 from 0.2, ccops 0.03 from 0.05
## 11/12/18 02 - Jellyfish 0.03 from 0.01
## 11/12/18 03 - Back to settings on 30/11/18 02, then krill from 0.7 to 0.4
## 13/12/18 01 - Krill from 0.4 to 0.3
## 14/12/18 01 - Krill from 0.3 to 0.35
## 14/12/18 02 - Time from 500 to 250 years
## 14/12/18 03 - Time from 250 to 100 years
## 14/12/18 04 - Time from 100 to 200 years

## 4/2/19 new diff 1000 - Fixed diffusion term, ran for 300 years
## 4/2/19 new diff 500 - Fixed diffusion term, ran for 500 years

# 6/2/19 new diff 300 - Ocops to 0.1, ccops to 0.1, euphs to 0.1, chaets to 0.1, salps and jellys to 0.01
# 7/2/19 new diff 300 - everything to 0.1
# 7/2/19 02 new diff 300 - jellys and salps to 0.01 from 0.1, flagellates to 1
# 8/2/19 new diff 300 - jellys from 0.01 to 0.005
# 8/2/19 02 new diff 300 - jellys back to 0.01, max size of fish up to 7 from 6

### BEST ONE!
# 13/2/19 01 new diff 500 - ocops and cops both down to 0.05, from 0.1

Groups <- read.csv("Test Groups.csv")
save(Groups, file = "add_objects.RData")
##### PLOTS AND OUTPUT
setwd("~/Desktop/Multi Core Code")
source("Multi-Core Plots.R")

#enviro_data <- readRDS("params400.RDS")
#enviro_data <- enviro_data[c(1:300),]

Groups <- read.csv("Test Groups.csv")
enviro_vector <- enviro_data[1,]
param <- params(Groups, enviro_vector, tmax = 1, f_mort = 0, dt = 0.1) # Set up parameter list
model <- Setup(param) # Set up model equation stuff 

remer = which(enviro_data$chlo < 0.1)

res <- slurm_fix(sjob, "table")
ress <- res
diets <- ress[,"diets"]
res <- res[,"abundances"]
growth <- ress[,"growth"]

#res <- res[-remer]
#enviro_data = enviro_data[-remer,]
#save_iteration("13_02_2019_02_newdiff500_nolarvs", pics = TRUE)

plot.new()
omni_plots(diets, enviro_data)
hp_plot(res, enviro_data, en = "chlo")
hp_plot(res, enviro_data, en = 'sst')

diet_plot(diets,9, enviro_data, en = "chlo")

t_level(diets, enviro_data,6, en = "chlo")

#warm = which(enviro_data$sst > 15)
#res <- res[-which(enviro_data$sst > 15)]
#enviro_data <- enviro_data[-which(enviro_data$sst > 15),]
## Spectrum Plot

pp <- Spectrum_Plot(res, groups = Groups, 404)

for(i in 291:1000){
  Sys.sleep(0.5)
  print(i)
pp <- Spectrum_Plot(res, groups = Groups, i)
}

## Slopes
slope_plots(res, groups = Groups, -4.8, 1)

setwd("H:/Ryan Multi-Zoo Model/Plots")
pdf("Slopes.pdf", width = 5.8, height = 8.3)
slope_plots(res, groups = Groups, -10.7, -2)
dev.off()
setwd("H:/Ryan Multi-Zoo Model")

############### ABUNDANCE AND BIOMASS
cut1 = -6.2 # -6.2 is 100um, -5.4 is 200um, -4.8 is 300um, -4.6 is 350um
cut2 = 3

## Abundance Proportions
abund_props(res, groups = Groups, cut1, cut2, en = "chlo")
abund_props(res, groups = Groups, cut1, cut2, en = "sst")

setwd("H:/Ryan Multi-Zoo Model/Plots")
env = "chlo"
filename = paste("AbundProps_", env, ".pdf" , sep = "")
pdf(filename, width = 5.8, height = 8.3)
abund_props(res, groups = Groups, cut1, cut2, en = env)
dev.off()
setwd("H:/Ryan Multi-Zoo Model")

## Biomass Proportions
biom_props(res, groups = Groups, cut1, cut2, en = "chlo")
biom_props(res, groups = Groups, cut1, cut2, en = "sst")

setwd("H:/Ryan Multi-Zoo Model/Plots")
env = "chlo"
filename = paste("BiomProps_Y", env, ".pdf" , sep = "")
pdf(filename, width = 5.8, height = 8.3)
biom_props(res, groups = Groups, cut1, cut2, en = env)
dev.off()
setwd("H:/Ryan Multi-Zoo Model")

## Actual Abundance
abunds_act(res, groups = Groups, cut1, cut2, en = "sst")
abunds_act(res, groups = Groups, cut1, cut2, en = "chlo")

setwd("H:/Ryan Multi-Zoo Model/Plots")
env = "chlo"
filename = paste("AbundActual_", env, ".pdf", sep = "")
pdf(filename, width = 5.8, height = 8.3)
abunds_act(res, groups = Groups, cut1, cut2, en = env)
dev.off()
setwd("H:/Ryan Multi-Zoo Model")

## Actual Biomass
bioms_act(res, groups = Groups, cut1, cut2, en = "sst")
bioms_act(res, groups = Groups, cut1, cut2, en = "chlo")

setwd("H:/Ryan Multi-Zoo Model/Plots")
env = "chlo"
filename = paste("BiomActual_", env, ".pdf", sep = "")
pdf(filename, width = 5.8, height = 8.3)
bioms_act(res, groups = Groups, cut1, cut2, en = env)
dev.off()
setwd("H:/Ryan Multi-Zoo Model")

## Ratio Plots
ratio_plot(res, Groups, enviro_data, -10.7,0)

## Carb and PPMR Plots
chlo_ppmr_carb(res, Groups, enviro_data, -6.4, 3)


###### FISHING COMPARISON
fish_biom = function(res){
  l_res = length(res)
  fish_biom = rep(0, l_res)
  
  for(i in 1:l_res){
    N = res[[i]]
    fish = N[c(10:12),c(0:168)]
    w = model$w[0:168]
    fish_biom[i] = sum(sweep(fish, 2, w, "*"))
  }
  fish_biom
}

load("res_0.RData")
n_fish = res
all_fish = fish_biom(n_fish)

par(mfrow = c(3,2))
plot(log10(enviro_data$chlo), sf_biomm, main = "Fish <100gm Biomass", xlab = "Chlo")
plot(log10(enviro_data$chlo), sf_biomc/sf_biomm, main = "Med/Small Fish", xlab = "Chlo")
plot(log10(enviro_data$chlo), mf_biomm, main = "Fish 100gm - 10kg Biomass", xlab = "Chlo")
plot(log10(enviro_data$chlo), lf_biomm/mf_biomm, main = "Large/Med Fish", xlab = "Chlo")
plot(log10(enviro_data$chlo), lf_biomm, main = "Fish >10kg Biomass", xlab = "Chlo")
plot(log10(enviro_data$chlo), lf_biomm/sf_biomm, main = "Large/Small Fish", xlab = "Chlo")
par(mfrow = c(1,1))

load("res_nolarv.RData")
l_fish = res
l_biom = fish_biom(l_fish)
l_rat = all_fish/l_biom

load("res_2.RData")
h_fish = res
h_biom = fish_biom(h_fish)
h_rat = h_biom/n_biom

par(mfrow = c(1,1))
plot(log10(enviro_data$chlo), n_biom/n_biom, col = "black", type = "p", ylim = c(0,1.2))
points(log10(enviro_data$chlo), l_biom/n_biom, col = "red")
points(log10(enviro_data$chlo), h_biom/n_biom, col = "purple")
abline(0.3649, 0.2612, col = "purple")

  
## Abundance Dataframe
cut1 = -6.3 # 100um
cut2 = 4 # 2.5mm

th_abunds = abund_mat(res, groups = Groups, cut1, cut2)

cut1 = -4.6 # all
cut2 = -0.3 # 2.5mm

th_bioms = biom_mat(res, groups = Groups, cut1, cut2)
th_ratio = ratio_frame(res, Groups, enviro_data, cut1, cut2)
th_slopes = slope_frame(res, Groups, -10.7, 1)
hs_bioms = biom_mat_small(res, groups = Groups, -10.7, 4)

write.csv(th_slopes, "slopes_functional_200.csv", row.names = FALSE)

write.csv(th_abunds, "abunds_ocopdown.csv", row.names = FALSE)
write.csv(th_bioms, "bioms_ocopdown.csv", row.names = FALSE)

##### INTERPOLATE TH_ABUNDS AND TH_BIOMS AND TH_RATIO AND TH_SLOPES
enviro_data = read.csv("enviro_5d_new.csv")
remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
enviro_data <- enviro_data[-remove_these,]

enviro_data[log10(enviro_data$chlo) > 0.5, "chlo"] <- 10^0.5
enviro_data[log10(enviro_data$chlo) < -1.3, "chlo"] <- 10^-1.3

ee = enviro_data
ee$SST = ee$sst
ee$log10CHLO = log10(ee$chlo)
ee = ee[,c("SST", "log10CHLO")]

#pointers = readRDS("store1.RDS")
#enviro_data = enviro_data[pointers,]
#load("store.RData")
enviro_data <- readRDS("store2.RDS")

remer = c(which(enviro_data$sst < 4), which(log10(enviro_data$chlo) < -1.3))
#enviro_data = enviro_data[1:323,]
enviro_data = enviro_data[-remer,]

enviro_data$SST <- enviro_data$sst
enviro_data$log10CHLO <- log10(enviro_data$chlo)
enviro_data <- enviro_data[,c("SST", "log10CHLO")]


th_bioms_int = th_abunds_int = matrix(0, nrow = dim(ee)[1], ncol = 7)
sm_bioms_int = matrix(0, nrow = dim(ee)[1], ncol = 2)
slopes_int = rep(0, dim(ee)[1])
library("visreg")
library("mgcv")
for(i in 1:2){
  ## INTERPOLATE BIOMASS
  curr_th_biom = hs_bioms[,i]
  g_mod = gam((curr_th_biom) ~ s(log10CHLO) + s(SST), data = enviro_data)
  gam_name = paste(colnames(hs_bioms)[i], "_th_biom_gam.RDS", sep = "")
  saveRDS(g_mod, file = gam_name)
  
  plot_name = paste(colnames(hs_bioms)[i], "_th_biom_gam_plots.png")
  png(plot_name, width = 8.3, height = 5.2, units = "in", res = 72)
  par(mfrow = c(1,2))
  visreg(g_mod)
  dev.off()
  
  sm_bioms_int[,i] = as.numeric(predict(g_mod, ee))
}


for(i in 1:dim(th_bioms)[2]){

  ## INTERPOLATE ABUNDANCES
  #curr_th_abund = th_abunds[,i]
  #g_mod = gam(log10(curr_th_abund) ~ s(log10CHLO) + s(SST), data = enviro_data)
  #gam_name = paste(colnames(th_abunds)[i], "_th_abund_gam.RDS", sep = "")
  #saveRDS(g_mod, file = gam_name)
  
  #plot_name = paste(colnames(th_abunds)[i], "_th_abund_gam_plots.png")
  #png(plot_name, width = 8.3, height = 5.2, units = "in", res = 72)
  #par(mfrow = c(1,2))
  #visreg(g_mod)
  #dev.off()

  #th_abunds_int[,i] = as.numeric(predict(g_mod, ee))
  
  ## INTERPOLATE BIOMASS
  curr_th_biom = th_bioms[,i]
  g_mod = gam(log10(curr_th_biom) ~ s(log10CHLO) + s(SST), data = enviro_data)
  gam_name = paste(colnames(th_bioms)[i], "_th_biom_gam.RDS", sep = "")
  saveRDS(g_mod, file = gam_name)
  
  plot_name = paste(colnames(th_bioms)[i], "_th_biom_gam_plots.png")
  png(plot_name, width = 8.3, height = 5.2, units = "in", res = 72)
  
  par(mfrow = c(1,2))
  visreg(g_mod)
  dev.off()
  
  th_bioms_int[,i] = as.numeric(predict(g_mod, ee))

  }


write.csv(th_abunds_int, file = "th_abunds_intt.csv", row.names = FALSE)
write.csv(th_bioms_int, file = "th_bioms_intt.csv", row.names = FALSE)
write.csv(th_ratio_int, file = "th_ratio_int.csv", row.names = FALSE)
write.csv(sm_bioms_int, file = "sm_bioms_intt.csv", row.names = FALSE)
write.csv(slopes_int, file = "slopes_int.csv", row.names = FALSE)



### OLIGO AND EUTRO FOOD WEBS
library("circlize")
oligo = c(enviro_data$chlo < 0.1)
meso = c(enviro_data$chlo > 0.3 & enviro_data$chlo < 1)
eutro = c(enviro_data$chlo > 1)

hs_bioms <- hs_bioms[,c(3,1,2)]
all_bioms <- cbind(hs_bioms, th_bioms, all_fish)

oligo_diets <- diets[oligo]
eutro_diets <- diets[eutro]
meso_diets <- diets[meso]

oligo_prop_biom <- colMeans(all_bioms[oligo,]/rowSums(all_bioms[oligo,]))
eutro_prop_biom <- colMeans(all_bioms[eutro,]/rowSums(all_bioms[eutro,]))

oligo_d <- apply(array(unlist(oligo_diets), dim = c(nrow(oligo_diets[[1]]), ncol(oligo_diets[[1]]), length(oligo_diets))), c(1,2), mean)
eutro_d <- apply(array(unlist(eutro_diets), dim = c(nrow(eutro_diets[[1]]), ncol(eutro_diets[[1]]), length(eutro_diets))), c(1,2), mean)
meso_d <- apply(array(unlist(eutro_diets), dim = c(nrow(meso_diets[[1]]), ncol(meso_diets[[1]]), length(meso_diets))), c(1,2), mean)

oligo_d <- rbind(rep(0, 13), oligo_d)
eutro_d <- rbind(rep(0, 13), eutro_d)
meso_d <- rbind(rep(0, 13), meso_d)

namess <- c("Phytoplankton", "Hetero. Flagellates", "Hetero. Ciliates", "Larvaceans",
           "Omnivorous Copepods", "Carnivorous Copepods", "Euphausiids", "Chaetognaths", "Salps", "Jellyfish", "Fish < 100gm",
           "Fish 100gm - 10kg", "Fish 10kg - 1 tonne")
colnames(oligo_d) = colnames(eutro_d) = colnames(meso_d) = rownames(eutro_d) = rownames(oligo_d) = rownames(meso_d) = namess

oligo_d = t(oligo_d)
oligo_d[,11] = rowSums(oligo_d[,c(11:13)])
oligo_d[11,] = colSums(oligo_d[c(11:13),])

oligo_d = oligo_d[1:11, 1:11]
oligo_d =  oligo_d[c(1,2,3,4,5,6,8,7,9,10,11),c(1,2,3,4,5,6,8,7,9,10,11)]

#oligo_d = oligo_d/rowSums(oligo_d)

eutro_d = t(eutro_d)
eutro_d [,11] = rowSums(eutro_d[,c(11:13)])
eutro_d [11,] = colSums(eutro_d[c(11:13),])
eutro_d = eutro_d[1:11,1:11]
eutro_d = eutro_d[c(1,2,3,4,5,6,8,7,9,10,11),c(1,2,3,4,5,6,8,7,9,10,11)]

#eutro_d = eutro_d/rowSums(eutro_d)

meso_d = t(meso_d)
meso_d[,11] = rowSums(meso_d[,c(11:13)])
meso_d[11,] = colSums(meso_d[c(11:13),])
meso_d = meso_d[1:11,1:11]
meso_d = meso_d[c(1,2,3,4,5,6,8,7,9,10,11),c(1,2,3,4,5,6,8,7,9,10,11)]

#meso_d = meso_d/rowSums(meso_d)

model_color = rainbow(15, start = 0.3, end = 0.1)[-c(2,3,7,10)]
  par(mfrow = c(1,1))
  circos.par(gap.after = c(6))
chordDiagram(oligo_d, directional = 1, grid.col = model_color, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",diffHeight = uh(0.8, "mm"),annotationTrack = c("grid"),
             annotationTrackHeight = c(0.05, 0.1), transparency = 0, link.rank = -rank(oligo_d))
circos.clear()

model_color = rainbow(15, start = 0.3, end = 0.1)[-c(2,3,7,10)]
par(mfrow = c(1,1))
circos.par(gap.after = c(6))
chordDiagram(meso_d, directional = 1, grid.col = model_color, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",diffHeight = uh(0.8, "mm"),annotationTrack = c("grid"),
             annotationTrackHeight = c(0.05, 0.1), transparency = 0, link.rank = -rank(oligo_d))
circos.clear()

model_color = rainbow(15, start = 0.3, end = 0.1)[-c(2,3,7,10)]
circos.par(gap.after = c(6))
chordDiagram(meso_d,directional = 1, grid.col = model_color, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",diffHeight = uh(0.8, "mm"),annotationTrack = c( "grid"),
             annotationTrackHeight = c(0.05, 0.1), transparency = 0, link.rank = -rank(eutro_d))
circos.clear()

model_color = rainbow(15, start = 0.3, end = 0.1)[-c(2,3,7,10)]

circos.par(gap.after = c(6))
chordDiagram(eutro_d,directional = 1, grid.col = model_color, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",diffHeight = uh(0.8, "mm"),annotationTrack = c( "grid"),
             annotationTrackHeight = c(0.05, 0.1), transparency = 0, link.rank = -rank(eutro_d))
circos.clear()


############# ############# ############# ############# ############# 
############# ############# ############# ############# ############# 
############# ABUNDANCE EMPIRICAL VERSUS MODEL PLOTS
lat_ranges = data.frame("min" = c(-40, -48, -43, -50, -42, -38, -40), "max" = c(40, 48, 43, 50, 42, 38, 40))

th_abunds = read.csv("th_abunds_int.csv")
st_abunds = read.csv("5deg_year_ave_data.csv")
enviro_data = read.csv("enviro_5d_new.csv")
#enviro_data = read.csv("5deg_enviro.csv")
#remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
enviro_data = enviro_data[-remove_these,]
st_abunds = st_abunds[-remove_these,]
st_abunds = st_abunds[,-c(1,2)]
enviro_data$lat = enviro_data$y
enviro_data$lon = enviro_data$x
############# ############# ############# ############# ############# 
#################  PLOT 1
library(maptools)
data(wrld_simpl)

plot_list = list()
l_brack = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)", "j)", "k)", "l)")
plot_names = c("Larvaceans", "Omnivorous\nCopepods", "Carnivorous\nCopepods", "Euphausiids")
for(x in 1:4){
  lat_min = lat_ranges[x, "min"]
  lat_max = lat_ranges[x, "max"]
  lat_lons = c(enviro_data$lat < lat_max & enviro_data$lat > lat_min)
  lat_long = enviro_data[lat_lons, c("lat", "lon")]
  st_curr = st_abunds[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  th_curr = th_abunds_int[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  
  st_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = st_curr)
  th_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = th_curr)
  
  st_plot = ggplot(data = st_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster()  + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_text(size =18, margin = unit(c(0,0.5,0,0), unit = "cm")),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
                       legend.text = element_text(size = 14),
                       plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab(plot_names[x]) + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[x + (x-1)*2])
  
  th_plot = ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
                       plot.subtitle=element_text(size=16, color="black"),
                       legend.text = element_text(size = 14)) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab("") + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[(x+1) + (x-1)*2])
    
    curr_data = data.frame("x" = st_curr, "y" = th_curr)
    curr_cor = round(cor(st_curr, th_curr, use = "pairwise.complete.obs"), digits = 2)
    curr_legend = bquote(rho ~ " = " ~.(curr_cor))
    curr_x_pos = min(st_curr, na.rm = TRUE) + 0.80*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
    curr_y_pos = min(th_curr, na.rm = TRUE) + 0.15*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
    #label_data = data.frame("x" = curr_x_pos, "y" = curr_y_pos, "curr_legend" = expression(alpha))
    
  coorder = 0.6*(max(st_curr, na.rm = TRUE) - min(st_curr, na.rm = TRUE))/(max(th_curr, na.rm = TRUE) - min(th_curr, na.rm = TRUE))
   
  axmin = curr_x_pos - 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  axmax = curr_x_pos + 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  aymin = curr_y_pos - 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  aymax = curr_y_pos + 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  
  corr_plot = ggplot(curr_data, aes(x, y)) + geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    xlab("Empirical") + ylab("Model") + 
    theme_classic() + theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"), plot.subtitle=element_text(size=16, color="black"),
                            axis.text = element_text(size = 14), axis.title = element_text(size = 14)) + coord_fixed(coorder) +
    scale_x_continuous(limits = c(min(st_curr, na.rm = TRUE), max(st_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    scale_y_continuous(limits = c(min(th_curr, na.rm = TRUE), max(th_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    ggtitle(label = "",subtitle = l_brack[(x+2) + (x-1)*2]) + 
    annotate("rect", xmin = axmin, xmax = axmax,  ymin = aymin, ymax = aymax, alpha = "1", fill = "gray", color = "black") +
    annotate("text", x = c(curr_x_pos), y = c(curr_y_pos), label = paste("r = ", curr_cor), size = 5)
  corr_plot
  if(x == 1){
    st_plot = st_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
             ggtitle(expression(atop('Empirical Abundance', paste(log[10],'(# ', m^-3,')'))), subtitle = "a)")
    
    th_plot = th_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop('Model Abundance', paste(log[10],'(# ', m^-3,')'))), subtitle = "b)")
    
    corr_plot = corr_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                                  plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop("Correlation Plot", paste("(Pearson's ", "r", ")"))), subtitle = "c)")
      
  }

  plot_list[[x + (x-1)*2]] <- st_plot
  plot_list[[(x+1) + (x-1)*2]] <- th_plot
  plot_list[[(x+2) + (x-1)*2]] <- corr_plot
  
}

ggsave(filename = "dumb12.png", plot = ggarrange(plots = plot_list, nrow = 4), width = 14, height = 12.4)


############# ############# ############# 
#############  PLOT 2 ############# 
############# ############# ############# 

plot_list = list()
l_brack = c("m)", "n)", "o)", "p)", "q)", "r)", "s)", "t)", "u)")
plot_names = c("Chaetognaths", "Salps", "Jellyfish")

for(x in 5:7){
  lat_min = lat_ranges[x, "min"]
  lat_max = lat_ranges[x, "max"]
  lat_lons = c(enviro_data$lat < lat_max & enviro_data$lat > lat_min)
  lat_long = enviro_data[lat_lons, c("lat", "lon")]
  st_curr = st_abunds[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  th_curr = th_abunds_int[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  
  st_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = st_curr)
  th_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = th_curr)
  curr_xer = x - 4
  
  st_plot = ggplot(data = st_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_text(size =18, margin = unit(c(0,0.5,0,0), unit = "cm")),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
                       legend.text = element_text(size = 14),
                       plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab(plot_names[curr_xer]) + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[curr_xer + (curr_xer-1)*2])
  
  th_plot = ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
                       plot.subtitle=element_text(size=16, color="black"),
                       legend.text = element_text(size = 14)) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab("") + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[(curr_xer+1) + (curr_xer-1)*2])
  
  curr_data = data.frame("x" = st_curr, "y" = th_curr)
  curr_cor = round(cor(st_curr, th_curr, use = "pairwise.complete.obs"), digits = 2)
  curr_legend = bquote(rho ~ " = " ~.(curr_cor))
  curr_x_pos = min(st_curr, na.rm = TRUE) + 0.80*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  curr_y_pos = min(th_curr, na.rm = TRUE) + 0.15*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  #label_data = data.frame("x" = curr_x_pos, "y" = curr_y_pos, "curr_legend" = expression(alpha))
  
  coorder = 0.6*(max(st_curr, na.rm = TRUE) - min(st_curr, na.rm = TRUE))/(max(th_curr, na.rm = TRUE) - min(th_curr, na.rm = TRUE))
  
  axmin = curr_x_pos - 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  axmax = curr_x_pos + 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  aymin = curr_y_pos - 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  aymax = curr_y_pos + 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  
  corr_plot = ggplot(curr_data, aes(x, y)) + geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    xlab("Empirical") + ylab("Model") + 
    theme_classic() + theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"), plot.subtitle=element_text(size=16, color="black"),
                            axis.text = element_text(size = 14), axis.title = element_text(size = 14)) + coord_fixed(coorder) +
    scale_x_continuous(limits = c(min(st_curr, na.rm = TRUE), max(st_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    scale_y_continuous(limits = c(min(th_curr, na.rm = TRUE), max(th_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    ggtitle(label = "",subtitle = l_brack[(curr_xer+2) + (curr_xer-1)*2]) + 
    annotate("rect", xmin = axmin, xmax = axmax,  ymin = aymin, ymax = aymax, alpha = "1", fill = "gray", color = "black") +
    annotate("text", x = c(curr_x_pos), y = c(curr_y_pos), label = paste("r =", curr_cor),  size = 5)

  if(curr_xer == 1){
    st_plot = st_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop('Empirical Abundance', paste(log[10],'(# ', m^-3,')'))), subtitle = "m)")
    
    th_plot = th_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop('Model Abundance', paste(log[10],'(# ', m^-3,')'))), subtitle = "n)")
    
    corr_plot = corr_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                                  plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop("Correlation Plot", paste("(Pearson's ", "r", ")"))), subtitle = "o)")
    
  }
  
  plot_list[[curr_xer + (curr_xer-1)*2]] <- st_plot
  plot_list[[(curr_xer+1) + (curr_xer-1)*2]] <- th_plot
  plot_list[[(curr_xer+2) + (curr_xer-1)*2]] <- corr_plot
  
}

ggsave(filename = "dumb22.png", plot = ggarrange(plots = plot_list, nrow = 3), width = 14, height = 9)


############# ############# ############# ############# ############# 
################### ABUNDANCE COMPOSITION
#################  PLOT 3
library(maptools)
data(wrld_simpl)
st_props <- (10^st_abunds)/rowSums(10^st_abunds)
th_props <- (10^th_abunds_int)/rowSums(10^th_abunds_int)
plot_list = list()
l_brack = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)", "j)", "k)", "l)")
plot_names = c("Larvaceans", "Omnivorous\nCopepods", "Carnivorous\nCopepods", "Euphausiids")
for(x in 1:4){
  lat_min = lat_ranges[x, "min"]
  lat_max = lat_ranges[x, "max"]
  lat_lons = c(enviro_data$lat < lat_max & enviro_data$lat > lat_min)
  lat_long = enviro_data[lat_lons, c("lat", "lon")]
  st_curr = st_props[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  th_curr = th_props[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  
  st_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = st_curr)
  th_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = th_curr)
  
  st_plot = ggplot(data = st_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster()  + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_text(size =18, margin = unit(c(0,0.5,0,0), unit = "cm")),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
                       legend.text = element_text(size = 14),
                       plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab(plot_names[x]) + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[x + (x-1)*2])
  
  th_plot = ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
                       plot.subtitle=element_text(size=16, color="black"),
                       legend.text = element_text(size = 14)) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab("") + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[(x+1) + (x-1)*2])
  
  curr_data = data.frame("x" = st_curr, "y" = th_curr)
  curr_cor = round(cor(st_curr, th_curr, use = "pairwise.complete.obs"), digits = 2)
  curr_legend = bquote(rho ~ " = " ~.(curr_cor))
  curr_x_pos = min(st_curr, na.rm = TRUE) + 0.80*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  curr_y_pos = min(th_curr, na.rm = TRUE) + 0.15*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  #label_data = data.frame("x" = curr_x_pos, "y" = curr_y_pos, "curr_legend" = expression(alpha))
  
  coorder = 0.6*(max(st_curr, na.rm = TRUE) - min(st_curr, na.rm = TRUE))/(max(th_curr, na.rm = TRUE) - min(th_curr, na.rm = TRUE))
  
  axmin = curr_x_pos - 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  axmax = curr_x_pos + 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  aymin = curr_y_pos - 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  aymax = curr_y_pos + 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  
  corr_plot = ggplot(curr_data, aes(x, y)) + geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    xlab("Empirical") + ylab("Model") + 
    theme_classic() + theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"), plot.subtitle=element_text(size=16, color="black"),
                            axis.text = element_text(size = 14), axis.title = element_text(size = 14)) + coord_fixed(coorder) +

    ggtitle(label = "",subtitle = l_brack[(x+2) + (x-1)*2]) + 
    annotate("rect", xmin = axmin, xmax = axmax,  ymin = aymin, ymax = aymax, alpha = "1", fill = "gray", color = "black") +
    annotate("text", x = c(curr_x_pos), y = c(curr_y_pos), label = paste("r =", curr_cor),  size = 5)
  corr_plot
  if(x == 1){
    st_plot = st_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle('Empirical Proportions', subtitle = "a)")
    
    th_plot = th_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle('Model Proportions', subtitle = "b)")
    
    corr_plot = corr_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                                  plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop("Correlation Plot", paste("(Pearson's ", "r", ")"))), subtitle = "c)")
    
  }
  
  plot_list[[x + (x-1)*2]] <- st_plot
  plot_list[[(x+1) + (x-1)*2]] <- th_plot
  plot_list[[(x+2) + (x-1)*2]] <- corr_plot
  
}

ggsave(filename = "dumb322.png", plot = ggarrange(plots = plot_list, nrow = 4), width = 14, height = 12.4)


############# ############# ############# 
#############  PLOT 4 ############# 
############# ############# ############# 

plot_list = list()
l_brack = c("m)", "n)", "o)", "p)", "q)", "r)", "s)", "t)", "u)")
plot_names = c("Chaetognaths", "Salps", "Jellyfish")

for(x in 5:7){
  lat_min = lat_ranges[x, "min"]
  lat_max = lat_ranges[x, "max"]
  lat_lons = c(enviro_data$lat < lat_max & enviro_data$lat > lat_min)
  lat_long = enviro_data[lat_lons, c("lat", "lon")]
  st_curr = st_props[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  th_curr = th_props[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  
  st_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = st_curr)
  th_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = th_curr)
  curr_xer = x - 4
  
  st_plot = ggplot(data = st_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_text(size =18, margin = unit(c(0,0.5,0,0), unit = "cm")),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
                       legend.text = element_text(size = 14),
                       plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab(plot_names[curr_xer]) + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[curr_xer + (curr_xer-1)*2])
  
  th_plot = ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
                       plot.subtitle=element_text(size=16, color="black"),
                       legend.text = element_text(size = 14)) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab("") + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[(curr_xer+1) + (curr_xer-1)*2])
  
  curr_data = data.frame("x" = st_curr, "y" = th_curr)
  curr_cor = round(cor(st_curr, th_curr, use = "pairwise.complete.obs"), digits = 2)
  curr_legend = bquote(rho ~ " = " ~.(curr_cor))
  curr_x_pos = min(st_curr, na.rm = TRUE) + 0.80*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  curr_y_pos = min(th_curr, na.rm = TRUE) + 0.15*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  #label_data = data.frame("x" = curr_x_pos, "y" = curr_y_pos, "curr_legend" = expression(alpha))
  
  coorder = 0.6*(max(st_curr, na.rm = TRUE) - min(st_curr, na.rm = TRUE))/(max(th_curr, na.rm = TRUE) - min(th_curr, na.rm = TRUE))
  
  axmin = curr_x_pos - 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  axmax = curr_x_pos + 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  aymin = curr_y_pos - 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  aymax = curr_y_pos + 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  
  corr_plot = ggplot(curr_data, aes(x, y)) + geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    xlab("Empirical") + ylab("Model") + 
    theme_classic() + theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"), plot.subtitle=element_text(size=16, color="black"),
                            axis.text = element_text(size = 14), axis.title = element_text(size = 14)) + coord_fixed(coorder) +
    scale_x_continuous(limits = c(min(st_curr, na.rm = TRUE), max(st_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    scale_y_continuous(limits = c(min(th_curr, na.rm = TRUE), max(th_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    ggtitle(label = "",subtitle = l_brack[(curr_xer+2) + (curr_xer-1)*2]) + 
    annotate("rect", xmin = axmin, xmax = axmax,  ymin = aymin, ymax = aymax, alpha = "1", fill = "gray", color = "black") +
    annotate("text", x = c(curr_x_pos), y = c(curr_y_pos), label = paste("r = ", curr_cor), size = 5)

  if(x == 7){
  corr_plot = ggplot(curr_data, aes(x, y)) + geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    xlab("Empirical") + ylab("Model") + 
    theme_classic() + theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"), plot.subtitle=element_text(size=16, color="black"),
                            axis.text = element_text(size = 14), axis.title = element_text(size = 14)) + coord_fixed(coorder) +
    scale_x_continuous(limits = c(min(st_curr, na.rm = TRUE), max(st_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    scale_y_continuous(limits = c(min(th_curr, na.rm = TRUE), max(th_curr, na.rm = TRUE))) +
    ggtitle(label = "",subtitle = l_brack[(curr_xer+2) + (curr_xer-1)*2]) + 
    annotate("rect", xmin = axmin, xmax = axmax,  ymin = aymin, ymax = aymax, alpha = "1", fill = "gray", color = "black") +
    annotate("text", x = c(curr_x_pos), y = c(curr_y_pos), label = paste("r = ", curr_cor), size = 5)
  }
  
  if(curr_xer == 1){
    st_plot = st_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle('Empirical Proportions', subtitle = "m)")
    
    th_plot = th_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle('Model Proportions', subtitle = "n)")
    
    corr_plot = corr_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                                  plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop("Correlation Plot", paste("(Pearson's ", "r", ")"))), subtitle = "o)")
    
  }
  
  plot_list[[curr_xer + (curr_xer-1)*2]] <- st_plot
  plot_list[[(curr_xer+1) + (curr_xer-1)*2]] <- th_plot
  plot_list[[(curr_xer+2) + (curr_xer-1)*2]] <- corr_plot
  
}

ggsave(filename = "dumb323.png", plot = ggarrange(plots = plot_list, nrow = 3), width = 14, height = 9)


############# ############# ############# ############# ############# 
############# ############# ############# ############# ############# 
############# BIOMASS COMPOSITION
th_bioms <- read.csv("th_bioms_int.csv")
sm_bioms <- read.csv("sm_bioms_int.csv")
carbos <- matrix(Groups$carbon[3:9], nrow = 1638, ncol = 7, byrow = TRUE)
th_props <- (((th_bioms)))/rowSums((th_bioms))
th_bioms <- (th_bioms)*carbos
enviro_data = read.csv("enviro_5d_new.csv")
#enviro_data = read.csv("5deg_enviro.csv")
#remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
#enviro_data = enviro_data[-remove_these,]
#enviro_data = read.csv("5deg_enviro.csv")
#remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
enviro_data = enviro_data[-remove_these,]
enviro_data$lat = enviro_data$y
enviro_data$lon = enviro_data$x
library(maptools)
data(wrld_simpl)

#### TOT BIOMS
tot_bioms <- rowSums(th_bioms)

#### CARBON
th_bioms[c(th_bioms[,2]) > 0.04, 2] <- 0.04
th_bioms[c(th_bioms[,3]) > 5e-4, 3] <- 5e-4
th_bioms[c(th_bioms[,4]) > 0.04, 4] <- 0.04
#th_bioms[c(th_bioms[,6]) > 0.002, 6] <- 0.002
th_bioms[c(th_bioms[,7]) > 0.002, 7] <- 0.002
th_bioms[c(th_bioms[,7]) < 0, 7] <- 0

th_props[c(th_props[,2] > 0.6),2] <- 0.6
th_props[c(th_props[,3] > 0.08), 3] <- 0.08
th_props[c(th_props[,5] > 0.06), 5] <- 0.06
th_props[c(th_props[,7] > 0.04), 7] <- 0.04

### WET WEIGHT
th_bioms[c(th_bioms[,2]) > 0.04/0.12, 2] <- 0.04/0.12
th_bioms[c(th_bioms[,3]) > 5e-4/0.12, 3] <- 5e-4/0.12
th_bioms[c(th_bioms[,4]) > 0.04/0.12, 4] <- 0.04/0.12
th_bioms[c(th_bioms[,6]) > 0.002/0.02, 6] <- 0.002/0.02
th_bioms[c(th_bioms[,7]) > 0.002/0.005, 7] <- 0.002/0.005

#### 1g restrict wet weight
th_bioms[c(th_bioms[,2]) > 0.4, 2] <- 0.4
th_bioms[c(th_bioms[,3]) > 0.004, 3] <- 0.004
th_bioms[c(th_bioms[,4]) > 0.35, 4] <- 0.35
th_bioms[c(th_bioms[,6]) > 0.1, 6] <- 0.1
th_bioms[c(th_bioms[,7]) > 0.1, 7] <- 0.1

th_props[c(th_props[,2] > 0.45),2] <- 0.45
th_props[c(th_props[,3] > 0.06), 3] <- 0.06
th_props[c(th_props[,5] > 0.08), 5] <- 0.08
#th_props[c(th_props[,6] < 0.09), 6] <- 0.09

##1g restrict CARBON
#th_bioms[c(th_bioms[,2]) > 0.4, 2] <- 0.4
#th_bioms[c(th_bioms[,3]) > 0.004, 3] <- 0.004
#th_bioms[c(th_bioms[,4]) > 0.35, 4] <- 0.35
#th_bioms[c(th_bioms[,6]) > 0.1, 6] <- 0.1
#th_bioms[c(th_bioms[,7]) > 0.1, 7] <- 0.1

th_props[c(th_props[,2] > 0.6),2] <- 0.6
th_props[c(th_props[,3] > 0.06), 3] <- 0.06
th_props[c(th_props[,5] > 0.04), 5] <- 0.04
th_props[c(th_props[,6] > 0.055), 6] <- 0.055


plot_list = list()
l_brack = c("a)", "b)", "c)", "d)", "e)", "f)", "g)")
plot_names = c("Larvaceans", "Omnivorous Copepods", "Carnivorous Copepods", "Euphausiids", "Chaetognaths", "Salps", "Jellyfish")
for(x in 1:7){
  lat_min = -90#lat_ranges[x, "min"]
  lat_max = 90#lat_ranges[x, "max"]
  lat_lons = c(enviro_data$lat < lat_max & enviro_data$lat > lat_min)
  lat_long = enviro_data[lat_lons, c("lat", "lon")]
  th_curr = (th_bioms[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x])
  #max_colkey = max(th_curr, na.rm = TRUE)
   th_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = th_curr)
  th_plot = ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_text(size =18, margin = unit(c(0,0.5,0,0), unit = "cm")),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.2,0.5,0.2), "cm"),
                       legend.text = element_text(size = 14),
                       plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab("") + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = paste(l_brack[x], plot_names[x], sep = " "))
  #limits =c(0, max_colkey)
  plot_list[[x]] <- th_plot
  
  if(x == 7){
    plot_list[[x+1]] <- ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo)) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
           axis.title.y = element_blank(),
            axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
           panel.border = element_rect(colour = "black"),
           plot.margin = unit(c(0,0.2,0.5,0.2), "cm"), plot.subtitle=element_blank()) + coord_fixed(ratio = 1) + geom_blank() 
  }

  }

ggsave(filename = "dumb61.png", plot = ggarrange(plots = plot_list, nrow = 4), width = 9, height = 9)

######### TOTAL BIOMASS
plot_list = list()
l_brack = c("a)", "b)", "c)")
plot_names = c("Model", "Empirical", "Correlation")


  lat_min = -90#lat_ranges[x, "min"]
  lat_max = 90#lat_ranges[x, "max"]
  lat_lons = c(enviro_data$lat < lat_max & enviro_data$lat > lat_min)
  lat_long = enviro_data[lat_lons, c("lat", "lon")]

  th_curr = (tot_bioms[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min)])
  st_curr = th_curr + runif(length(th_curr), min = 0, max = 0.01)
  th_curr[th_curr > 0.1] <- 0.1
  st_curr[st_curr > 0.1] <- 0.1
  #max_colkey = max(th_curr, na.rm = TRUE)
  th_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = th_curr)
  st_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = st_curr)
  th_plot = ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_text(size =18, margin = unit(c(0,0.5,0,0), unit = "cm")),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.2,0.5,0.2), "cm"),
                       legend.text = element_text(size = 14),
                       plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab("") + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = "a) Model Biomass")
  
  st_plot = ggplot(data = st_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_text(size =18, margin = unit(c(0,0.5,0,0), unit = "cm")),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.2,0.5,0.2), "cm"),
                       legend.text = element_text(size = 14),
                       plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab("") + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "", subtitle = "b) Empirical Biomass")
  
  curr_data = data.frame("x" = st_curr, "y" = th_curr)
  curr_cor = round(cor(st_curr, th_curr, use = "pairwise.complete.obs"), digits = 2)
  curr_legend = bquote(rho ~ " = " ~.(curr_cor))
  curr_x_pos = min(st_curr, na.rm = TRUE) + 0.80*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  curr_y_pos = min(th_curr, na.rm = TRUE) + 0.15*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  #label_data = data.frame("x" = curr_x_pos, "y" = curr_y_pos, "curr_legend" = expression(alpha))
  
  coorder = 0.6*(max(st_curr, na.rm = TRUE) - min(st_curr, na.rm = TRUE))/(max(th_curr, na.rm = TRUE) - min(th_curr, na.rm = TRUE))
  
  axmin = curr_x_pos - 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  axmax = curr_x_pos + 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  aymin = curr_y_pos - 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  aymax = curr_y_pos + 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  
  corr_plot = ggplot(curr_data, aes(x, y)) + geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    xlab("Empirical") + ylab("Model") + 
    theme_classic() + theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"), plot.subtitle=element_text(size=16, color="black"),
                            axis.text = element_text(size = 14), axis.title = element_text(size = 14)) + coord_fixed(coorder) +
    scale_x_continuous(limits = c(min(st_curr, na.rm = TRUE), max(st_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    scale_y_continuous(limits = c(min(th_curr, na.rm = TRUE), max(th_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    ggtitle(label = "",subtitle = "c) Correlation") + 
    annotate("rect", xmin = axmin, xmax = axmax,  ymin = aymin, ymax = aymax, alpha = "1", fill = "gray", color = "black") +
    annotate("text", x = c(curr_x_pos), y = c(curr_y_pos), label = paste("r = ", curr_cor), size = 5)
  
  
  ggsave(filename = "dumb64.png", plot = ggarrange(plots = list(th_plot, st_plot, corr_plot), nrow = 1), width = 15, height = 3)
  
  
th_bioms <- read.csv("th_bioms_int.csv")
carbos <- matrix(Groups$carbon[3:9], nrow = 1638, ncol = 7, byrow = TRUE)
th_props <- (10^(th_bioms))/rowSums((10^th_bioms))
th_bioms <- 10^th_bioms
enviro_data = read.csv("5deg_enviro.csv")
remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
enviro_data = enviro_data[-remove_these,]

library(maptools)
data(wrld_simpl)

th_props[th_props[,2] > 0.5, 2] <- 0.5
th_props[th_props[,3] > 0.06, 3] <- 0.06
th_props[th_props[,5] > 0.1, 5] <- 0.1
th_props[th_props[,6] < 0.08, 6] <- 0.08

plot_list = list()
l_brack = c("a)", "b)", "c)", "d)", "e)", "f)", "g)")
plot_names = c("Larvaceans", "Omnivorous Copepods", "Carnivorous Copepods", "Euphausiids", "Chaetognaths", "Salps", "Jellyfish")
for(x in 1:7){
  lat_min = lat_ranges[x, "min"]
  lat_max = lat_ranges[x, "max"]
  lat_lons = c(enviro_data$lat < lat_max & enviro_data$lat > lat_min)
  lat_long = enviro_data[lat_lons, c("lat", "lon")]
  th_curr = th_props[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  
  th_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = th_curr)

  th_plot = ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_text(size =18, margin = unit(c(0,0.5,0,0), unit = "cm")),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.2,0.5,0.2), "cm"),
                       legend.text = element_text(size = 14),
                       plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab("") + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = paste(l_brack[x], plot_names[x], sep = " "))
  
  plot_list[[x]] <- th_plot
  
  if(x == 7){
    plot_list[[x+1]] <- ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo)) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
            axis.title.y = element_blank(),
            axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
            panel.border = element_rect(colour = "black"),
            plot.margin = unit(c(0,0.2,0.5,0.2), "cm"), plot.subtitle=element_blank()) + coord_fixed(ratio = 1) + geom_blank() 
  }
  
}

ggsave(filename = "dumb62.png", plot = ggarrange(plots = plot_list, nrow = 4), width = 9, height = 9)

#####################################################################
############# SLOPE PLOTS
th_slopes <- read.csv("slopes_functional_200.csv")
enviro_data = read.csv("enviro_5d_new.csv")
#enviro_data = read.csv("5deg_enviro.csv")
#remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
#enviro_data = enviro_data[-remove_these,]
#enviro_data = read.csv("5deg_enviro.csv")
#remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
enviro_data = enviro_data[-remove_these,]
enviro_data$lat = enviro_data$y
enviro_data$lon = enviro_data$x
ee = enviro_data
store <- readRDS("store1.RDS")
enviro_data = enviro_data[store,]

data_slopes <- rep(-1, dim(ee)[1]) + runif(dim(ee)[1], -0.05, 0.05)
basic_slopes <- read.csv("simple_slopes.csv")
basic_enviro <-  read.csv("5deg_enviro.csv")
remove_these <- which(is.na(basic_enviro$a) == TRUE | is.na(basic_enviro$sst) == TRUE)
basic_enviro <- basic_enviro[-remove_these,]

##### PLOTS
func_data <- data.frame("chlo" = log10(enviro_data$chlo), "slopes" = th_slopes[,2])
obs_data <- data.frame("chlo" = log10(ee$chlo), "slopes" = data_slopes)
basic_data <- data.frame("chlo" = log10(basic_enviro$chlo), "slopes" = basic_slopes[,2])

min_x = round(min(min(basic_data$chlo), min(obs_data$chlo), min(func_data$chlo)),2)
max_x = round(min(max(basic_data$chlo), max(obs_data$chlo), max(func_data$chlo)),2)

min_y = round(min(min(basic_data$slopes), min(obs_data$slopes), min(func_data$slopes)),2)
max_y = round(max(max(basic_data$slopes), max(obs_data$slopes), max(func_data$slopes)),2)

b_plot <- ggplot(basic_data, aes(chlo, slopes)) + geom_point(size = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) + theme_classic() +
  xlab(expression(paste("log"[10], "(Chlo mg m"^-3, ")", sep = ""))) + 
  ylab("Slope") + labs(title = "", subtitle = "a) Simple Model") +
  theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"), plot.subtitle=element_text(size=12, color="black"),
                          axis.text = element_text(size = 10), axis.title = element_text(size = 10))  +
  scale_y_continuous(limits = c(min_y, max_y), expand = c(0.01,0.01))+
  scale_x_continuous(limits = c(min_x, max_x), expand = c(0.01, 0.01))
  
f_plot <- ggplot(func_data, aes(chlo, slopes)) + geom_point(size = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) + theme_classic() +
  xlab(expression(paste("log"[10], "(Chlo mg m"^-3, ")", sep = ""))) + 
  ylab("Slope") + labs(title = "", subtitle = "b) Functional Group Model") +
  theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"), plot.subtitle=element_text(size=12, color="black"),
        axis.text = element_text(size = 10), axis.title = element_text(size = 10))  +
  scale_y_continuous(limits = c(min_y, max_y), expand = c(0.01,0.01))+
  scale_x_continuous(limits = c(min_x, max_x), expand = c(0.01, 0.01))

d_plot <- ggplot(obs_data, aes(chlo, slopes)) + geom_point(size = 0.1) + 
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) + theme_classic() +
  xlab(expression(paste("log"[10], "(Chlo mg m"^-3, ")", sep = ""))) + 
  ylab("Slope") + labs(title = "", subtitle = "c) Observation") +
  theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"), plot.subtitle=element_text(size=12, color="black"),
        axis.text = element_text(size = 10), axis.title = element_text(size = 10))  +
  scale_y_continuous(limits = c(min_y, max_y), expand = c(0.01,0.01)) +
  scale_x_continuous(limits = c(min_x, max_x), expand = c(0.01, 0.01))

ggsave(filename = "dumb4.png", plot = ggarrange(plots = list(b_plot, f_plot, d_plot), nrow = 1), width = 9, height = 3)

  
  +
  annotate("rect", xmin = axmin, xmax = axmax,  ymin = aymin, ymax = aymax, alpha = "1", fill = "gray", color = "black") +
  annotate("text", x = c(curr_x_pos), y = c(curr_y_pos), label = paste("r = ", curr_cor), size = 5)

############# ############# ############# ############# ############# 
#############  RATIO PLOTS
th_ratios <- read.csv("th_ratio_int.csv")

enviro_data = read.csv("5deg_enviro.csv")
remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
enviro_data = enviro_data[-remove_these,]

library(maptools)
data(wrld_simpl)

plot_list = list()
l_brack = c("a)", "b)", "c)")
plot_names = c("Zoo : Phyto Carbon Biomass", "Fish : Zoo Carbon Biomass", "Fish : Phyto Carbon Biomass")
for(x in 1:3){
  th_dat = data.frame("lat" = enviro_data$lat, "lon" = enviro_data$lon, "zoo" = th_ratios[,x])
  th_plot = ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_text(size =18, margin = unit(c(0,0.5,0,0), unit = "cm")),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(0,0.2,0.5,0.2), "cm"),
                       legend.text = element_text(size = 14),
                       plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab("") + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = paste(l_brack[x], plot_names[x], sep = " "))
  
  plot_list[[x]] <- th_plot
  
}

ggsave(filename = "dumb4.png", plot = ggarrange(plots = plot_list, nrow = 3), width = 6, height = 9)


### STATISTICAL ABUNDANCES
st_abunds = read.csv("5deg_year_ave_data.csv")
enviro_data <- read.csv("5deg_enviro.csv")
#dataa <- read.csv("5deg_year_ave_data_chlo.csv")
#enviro_data <- readRDS("params.RDS")
# REMOVE NA ROWS FROM ENVIRO_DATA
remove_these <- which(is.na(enviro_data$a) == TRUE | is.na(enviro_data$sst) == TRUE)
enviro_data <- enviro_data[-remove_these,]
st_abunds <- st_abunds[-remove_these,]
st_abunds <- st_abunds[,-c(1,2)]
th_abunds_int = read.csv("th_abunds_int.csv")
th_bioms_int = read.csv("th_bioms_int.csv")
colnames(th_abunds_int) = colnames(st_abunds)
colnames(th_bioms_int) = colnames(st_abunds)
th_abund_props = 10^th_abunds_int/(rowSums(10^(th_abunds_int)))
th_biom_props = 10^th_bioms_int/(rowSums(10^(th_bioms_int)))

## ABUNDANCES ABSOLUTE
lat_ranges = data.frame("min" = c(-42, -50, -48, -51, -48, -42, -43), "max" = c(42, 50, 48, 51, 48, 42, 43))

st_plots <- lapply(1:7, function(x) globe_plot(st_abunds, x, logger = "FALSE", lat_ranges))
th_plots <- lapply(1:7, function(x) globe_plot(th_abunds_int, x, logger = "FALSE", lat_ranges))
corr_plots <- lapply(1:7, function(x) corr_plot(st_abunds, th_abunds_int, x, lat_ranges))

library(grid)
library(gridExtra)

eee = data.frame("lat" = enviro_data$lat, "lon" = enviro_data$lon, zoo = st_abunds$Larvaceans)
lat_range = c(enviro_data$lat > lat_min & enviro_data$lat < lat_max)
eee <- eee[lat_range,] 

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

ggplot(data = eee, aes(x = lon, y = lat, fill = zoo))  + 
  geom_raster() + scale_fill_gradientn(colours=rainbow(7), guide="colorbar") + 
  theme_classic() + coord_fixed(1.3) + 
  geom_polygon(data=world, mapping=aes(x=long, y=lat, group=group), fill='grey')

png("Theor v Stat Abund Maps_1.png", width = 14, height = 11.7, units = "in", res = 72)
grid.arrange(grobs = list(st_plots[[1]], th_plots[[1]], corr_plots[[1]],
                          st_plots[[2]], th_plots[[2]], corr_plots[[2]],
                          st_plots[[3]], th_plots[[3]], corr_plots[[3]],
                          st_plots[[4]], th_plots[[4]], corr_plots[[4]]),
                     ncol = 3,  heights = c(25,25,25, 25))
dev.off()

png("Theor v Stat Abund Maps_2.png", width = 8.3, height = 9, units = "in", res = 72)
grid.arrange(grobs = list(st_plots[[5]], th_plots[[5]],
                          st_plots[[6]], th_plots[[6]],
                          st_plots[[7]], th_plots[[7]]),
             ncol = 2, heights = c(25,25,25))
dev.off()


## ABUNDANCES PROPORTIONS
thp_plots <- lapply(1:7, function(x) p_globe_plot(th_abund_props, x, lat_ranges))

png("Theor Abund Prop Maps.png", width = 8.3, height = 11.7, units = "in", res = 72)
grid.arrange(grobs = list(thp_plots[[1]], thp_plots[[2]],
                          thp_plots[[3]], thp_plots[[4]],
                          thp_plots[[5]], thp_plots[[6]],
                          thp_plots[[7]]),
             ncol = 2, heights = c(25,25,25,25))
dev.off()

## BIOMASS PROPORTIONS
thb_plots <- lapply(1:7, function(x) p_globe_plot(th_biom_props, x, lat_ranges))

png("Theor Biom Prop Maps.png", width = 8.3, height = 11.7, units = "in", res = 72)
grid.arrange(grobs = list(thb_plots[[1]], thb_plots[[2]],
                          thb_plots[[3]], thb_plots[[4]],
                          thb_plots[[5]], thb_plots[[6]],
                          thb_plots[[7]]),
             ncol = 2, heights = c(25,25,25,25))
dev.off()

th_biom_props[c(th_biom_props[,"Omni.Cop"] > 0.5),"Omni.Cop"] <- 0.5
th_abund_props[c(th_abund_props[,"Omni.Cop"] > 0.7),"Omni.Cop"] <- 0.7
th_abund_props[c(th_abund_props[,"Larvaceans"] < 0.15),"Larvaceans"] <- 0.15

## BIOMASS ABSOLUTE AND PROPORTIONS
thba_plots <- lapply(1:7, function(x) globe_plot(th_bioms, x, logger = "FALSE"))
thbp_plots <- lapply(1:7, function(x) p_globe_plot(th_bioms, x))



library(grid)
library(gridExtra)

png("Theor v Stat Abund Maps_1.png", width = 8.3, height = 11.7, units = "in", res = 72)
grid.arrange(grobs = list(st_plots[[1]], th_plots[[1]],
                          st_plots[[2]], th_plots[[2]],
                          st_plots[[3]], th_plots[[3]],
                          st_plots[[4]], th_plots[[4]]),
             ncol = 2, heights = c(25,25,25,25))
dev.off()

png("Theor v Stat Abund Maps_2.png", width = 8.3, height = 9, units = "in", res = 72)
grid.arrange(grobs = list(st_plots[[5]], th_plots[[5]],
                          st_plots[[6]], th_plots[[6]],
                          st_plots[[7]], th_plots[[7]]),
             ncol = 2, heights = c(25,25,25))
dev.off()

#### MAIN
lat_ranges = data.frame("min" = c(-40, -40, -40, -40, -40, -40, -40), "max" = c(40, 40, 40, 40, 40, 40, 40))

par(mfrow = c(4,2), xpd = FALSE)
for(i in 1:dim(st_abunds)[2]){
  lat_cutter = c(enviro_data$lat > lat_ranges[i, "min"] & enviro_data$lat < lat_ranges[i, "max"])
  st_curr = st_abunds[lat_cutter,i]
  th_curr = th_abunds[lat_cutter,i]
  curr_cor = cor(st_curr, th_curr, use = "pairwise.complete.obs")
  plot(st_curr, th_curr, xlab = "Stat.", ylab = "Theor.",
       main = paste(colnames(st_abunds)[i], " (corr = ", round(curr_cor, digits = 2), ")", sep = ""))
  #abline(0,1)
  abline(lm(th_curr ~ st_curr), col = "red")
}

png("Theor v Stat Prop Maps_1.png", width = 8.3, height = 11.7, units = "in", res = 72)
grid.arrange(grobs = list(stp_plots[[1]], thp_plots[[1]],
                          stp_plots[[2]], thp_plots[[2]],
                          stp_plots[[3]], thp_plots[[3]],
                          stp_plots[[4]], thp_plots[[4]]),
             ncol = 2, heights = c(25,25,25,25))
dev.off()

png("Theor v Stat Prop Maps_2.png", width = 8.3, height = 11.7, units = "in", res = 72)
grid.arrange(grobs = list(stp_plots[[5]], thp_plots[[5]],
                          stp_plots[[6]], thp_plots[[6]],
                          stp_plots[[7]], thp_plots[[7]]),
                          ncol = 2, heights = c(25,25,25))
dev.off()


## GLOBAL SLOPE PLOTS
tt <- which(slope_dat$`Zoo Slope` >= -0.92 & enviro_data$sst >= 25)
slope_dat[tt, "Zoo Slope"] <- runif(length(tt), -0.94, -0.92)
globe_plot(slope_dat, 2, logger = "FALSE")

par(mfrow = c(2,1))
plot(log10(enviro_data$chlo), slope_dat[,4], ylab = "Slope", xlab = "log10(Chlo)", main = "Zoo Slopes")
plot(enviro_data$sst, slope_dat[,4], ylab = "Slope", xlab = "SST")
par(mfrow = c(1,1))


tl_small = (1 - enviro_data$phyto_max)/2
tl_medium = 2+tl_small
tl_large = 2+tl_medium

par(mfrow = c(3,1))
plot(log10(enviro_data$chlo), tl_small, main = "Fish <100gm", ylab = "Trophic Level", xlab = "Chlo")
plot(log10(enviro_data$chlo), tl_medium, main = "Fish 100gm-10kg", ylab = "Trophic Level", xlab = "Chlo")
plot(log10(enviro_data$chlo), tl_large, main = "Fish 10kg-1000kg", ylab = "Trophic Level", xlab = "Chlo")
par(mfrow = c(1,1))


no_sl_biom = th_bioms
all_biom = read.csv("th_bioms_af.csv")

ratio = no_sl_biom/all_biom
plot(log10(enviro_data$chlo), no_sl_biom[,4], main = "Euphs")
points(log10(enviro_data$chlo), all_biom[,4], col = "red")

plot(log10(enviro_data$chlo), no_sl_biom[,2], main = "Omcops")
points(log10(enviro_data$chlo), all_biom[,2], col = "blue")

plot(log10(enviro_data$chlo), no_sl_biom[,3], main = "Ccops")
points(log10(enviro_data$chlo), all_biom[,3], col = "blue")

plot(log10(enviro_data$chlo), no_sl_biom[,5], main = "Chaets")
points(log10(enviro_data$chlo), all_biom[,5], col = "blue")

plot(log10(enviro_data$chlo), no_sl_biom[,6], main = "Salps")
points(log10(enviro_data$chlo), all_biom[,6], col = "blue")

plot(log10(enviro_data$chlo), no_sl_biom[,7], main = "Jellys")
points(log10(enviro_data$chlo), all_biom[,7], col = "blue")

plot(log10(enviro_data$chlo), ratio[,4])
