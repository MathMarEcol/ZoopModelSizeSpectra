## An extension of the model published in Heneghan et al., (2016):
## Models multiple zooplankton functional groups, and a single fish community 

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
  
  cc_phyto = 0.15   # Carbon content of phytoplankton size classes
  
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
  assim_phyto =  (param$groups$alpha)*cc_phyto/param$nutrition # Phytoplankton
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
  diff_log_t_phyto = ((w^-1) %*% t(w_phyto^2))/log(10) # Diffusion
  diet_log_t_phyto = matrix(w_phyto, nrow = length(w), ncol = length(w_phyto), byrow = TRUE) # Diet/Ingestion
  
  # Predators are rows, dynam prey weights are columns
  gg_log_t_dynam = ((w^-1) %*% t(w))/log(10) # Growth
  diff_log_t_dynam = ((w^-1) %*% t(w^2))/log(10) # Diffusion
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
  temp_zoo <- rep(2.^((environ$sst - 15.6)/10), num_zoo) # exp(23.93 - 0.59/(8.62e-05*(273+environ$sst)))
  temp_fish <- rep(2.^((environ$sst - 15.6)/10), num_fish)
  temp_effect <- matrix(c(temp_zoo, temp_fish), nrow = ngrps, ncol = ngrid)
  
  #### CALCULATES CONSTANT BITS OF THE MODEL FUNCTIONS FOR EACH GROUP
  for(i in 1:ngrps){ 
    ## Senescence mortality
    model$M_sb[i,] = param$ZSpre*(w/(10^(grp$Wmat[i])))^param$ZSexp
    model$M_sb[i, 10^(grp$Wmax[i]) < w] = 0
    model$M_sb[i, 10^(grp$Wmat[i]) > w] = 0	
    
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
  model$M_sb[c(no_sen, ngrps)] = 0 
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
    A[,idx] <- dt/dx*(gg[,idx-1] + diff[,idx-1]/(2*dx))
    B[,idx] <- diff[,idx+1]*dt/(2*dx^2)
    C[,idx] <- 1 + dt*Z[,idx] + dt/dx*(gg[,idx] + diff[,idx]/dx)
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

Diet_Matrix = function(models, params){
  model = modelss
  param = params
  diet = model$diet
  fish_diet = model$fish_diet
  ngrps = dim(param$groups)[1]
  
  assim_dynam = matrix(param$groups$alpha*param$nutrition, nrow = ngrps, ncol = ngrps, byrow = TRUE)/
                    matrix(param$nutrition, nrow = ngrps, ncol = ngrps)
  
  tot_phyto_prod = ((10^1.58)*(param$environ$chlo^1.29)*1e-2*365)
  
  diet_matt = matrix(0, nrow = ngrps, ncol = (ngrps+1))
  
 #  fish_s = seq(1,which(round(log10(model$w), digits = 2) == 2),1) # (10^-3g - 10^2g)
 # fish_l = seq((length(fish_s)+1), length(model$w), 1) # (10^2g - 10^6g)
  
 # fish_diet_mat = colMeans(fish_diet[ceiling(0.5*dim(fish_diet)[1]):dim(fish_diet)[1],,], dim = 1)
  diet_mat = colMeans(diet[ceiling(0.5*dim(diet)[1]):dim(diet)[1],,,], dim = 1)
  diet_mat[1:ngrps, ,2:(ngrps+1)] = sweep(diet_mat[1:ngrps,,2:(ngrps+1)], c(1,3), assim_dynam, "/")
  diet_mat = aperm(diet_mat, c(2,1,3))
  
  diet_matt[(1:ngrps),1:(ngrps+1)] = colSums(diet_mat)

   ### DON'T NEED IF WE HAVE SMALL, MEDIUM, LARGE FISH
  # diet_matt[10,1:11] = colSums(diet_mat[fish_s,10,])
 #  diet_matt[11,1:11] = colSums(diet_mat[fish_l,10,])
 #  diet_matt[10,11] = fish_diet_mat[1,1]
 # diet_matt[11,11] = fish_diet_mat[2,1]
#  diet_matt[11,12] = fish_diet_mat[2,2]
    
#  diet_mat[11,,] = diet_mat[10,,]*matrix(fish_large, nrow = length(fish_large), ncol = 11)  

  diet_matt = round(diet_matt, digits = 6)*1000
  group_names = as.character(param$groups$species)
  rownames(diet_matt) = as.character(param$groups$species)
  colnames(diet_matt) = c("Phytoplankton", group_names)
  return(list(diet_matt, tot_phyto_prod))

}

#### PLOT AVERAGE SPECTRUM
Spectrum_Plot = function(models, params, fish_on){
  modelss = models
  param = params
  
  N_sav = modelss$N
  N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
  zoo_groups =  param$zoo_grps
  fish_groups = param$fish_grps
  num_zoo = param$num_zoo
  num_fish = param$num_fish
  
  if(length(zoo_groups) > 1){
  tot_zoo = colSums(N_ave[zoo_groups,])
  }
  
  if(length(zoo_groups) == 1){
  tot_zoo = N_ave[zoo_groups,]  
  }
  
  if(num_fish > 0){
  if(num_fish > 1){
    tot_fish = colSums(N_ave[fish_groups,])
  }else{tot_fish = N_ave[fish_groups,]}
  }
  
  y_d_ref = which(round(log10(model$w), digits = 2) == 5)
  if(fish_on == TRUE){
  y_down = floor(log10(N_ave[dim(N_ave)[1],][y_d_ref]))
  }else{y_down = floor(log10(N_ave[num_zoo,][which(model$w == 10)]))}
  y_up = ceiling(log10(model$nPP[5]))
  
  ## PLOT PHYTO-ZOO-FISH TOTALS
  par(mar = c(4,4,4,2))
  plot(log10(model$w_phyto), log10(model$nPP), type= "l", col = "green",
       lwd = 2, xlim = c(log10(model$w_phyto[5]), log10(model$w[length(model$w)])), 
       ylim = c(y_down, y_up), 
       xlab = "", 
       ylab ="")
  lines(log10(model$w), log10(tot_zoo), col = "red", lwd = 2)
  if(num_fish > 0){
  lines(log10(model$w), log10(tot_fish), col = "blue", lwd = 2)
  }
  mtext(expression(paste("log"[10], "(Abundance ", "# m "^-3, ")")), side = 2, line = 2.2, cex=1)
  mtext(expression(paste("log"[10], "(Body Weight, g)")), side = 1, line = 2.2, cex = 1)
  legend("bottomleft", legend = c("Phytoplankton", "Total Zooplankton", "Fish Community"),
         col = c("green","red", "blue"), lty = 1, lwd = 2, bty = "n", cex = 0.9)
  title(main = paste("Phytoplankton Slope = " , round(enviro$b, digits = 2),
                               "\nTemperature = ", round(enviro$sst), "C"))
  ## PLOT ZOO GROUPS
  if(num_zoo > 1){
  coll = rainbow(num_zoo)
  
  for(i in 1:num_zoo){
    colll = coll[i]
    lines(log10(model$w), log10(N_ave[zoo_groups[i],]), lty = 2, col = colll, lwd = 1.5)
  }
  legend(x = 0, y = y_up + 1.5, legend = as.character(param$groups$species[zoo_groups]),
         col = coll, lty = 2, lwd = 1.5,bty = "n", cex = 0.8)
  }
  ## PLOT FISH GROUPS
  if(num_fish > 0){
  for(i in 1:num_fish){
    lines(log10(model$w), log10(N_ave[fish_groups[i],]), lty = 2, col = "blue", lwd = 1.5)
  }
  }
  }

## PLOT BIOMASS OVER TIME
Biomass_Plot = function(models, params){
  modelss = models
  param = params
  
  zoo_groups =  param$zoo_grps
  fish_groups = param$fish_grps
  num_zoo = param$num_zoo
  num_fish = param$num_fish
  
  if(num_fish > 0){
  FB_total = modelss$Biomass[(ceiling(0.5*dim(modelss$Biomass)[1])):(dim(modelss$Biomass)[1]), fish_groups]
  }
  Z_groups = modelss$Biomass[(ceiling(0.5*dim(modelss$Biomass)[1])):(dim(modelss$Biomass)[1]), zoo_groups]
  
  if(num_zoo > 1){
  ZB_total = rowSums(Z_groups) 
  }else{
  ZB_total = Z_groups  
  }
  
  if(num_fish > 0){
  if(num_fish > 1){
    FB_total = rowSums(FB_total)
  }
  }
  
  par(mar = c(4,4,4,10))
  if(num_fish > 0){
  plot(FB_total, type = "l", col = "blue", lwd = 2, ylim = c(0, (max(c(FB_total, ZB_total)))), xaxt = "n", 
       ylab = "", xlab = "")
  }else{
    plot(0, type = "n", xlim = c(0, length(ZB_total)), ylim = c(0, max(ZB_total)), xaxt = "n", 
         ylab = "", xlab = "")
  }
  lines(ZB_total, type = "l", col = "red", lwd = 2)
  mtext(text = expression(paste("Total Biomass ", "g m"^-3)), side = 2, line = 2.5)
  mtext(text = "Time (years)", side = 1, line = 2.5)
  axis(side = 1, at = seq(1,length(ZB_total), length.out = 10), labels = round(seq(0.1*length(ZB_total)/12.17,length(ZB_total)/12.17,length.out =10)), cex.axis = 1)
  
  if(length(zoo_groups) > 1){
  coll = rainbow(num_zoo)
  for(i in 1:length(zoo_groups)){
    lines(Z_groups[,i], lty = 2, col = coll[i], lwd = 1)
  }
  }
  
  title(main = paste("Phytoplankton Slope = " , round(enviro$b, digits = 2),
                     "\nTemperature = ", round(enviro$sst), "C"))
  par(xpd = TRUE)
  max_x = length(ZB_total)
  
  if(length(zoo_groups) > 1){
  if(num_fish > 0){
  max_y = max(c(FB_total, ZB_total))
  legend(x= (max_x + 2), y = max_y, legend = c("All Zooplankton", "All Fish", as.character(param$groups$species[zoo_groups])),
         col = c("Red", "Blue", coll), lty = c(1,1, rep(2, length(zoo_groups))), lwd = c(2,2, rep(1.5, length(zoo_groups))),bty = "n", cex = 0.8)
  par(xpd = FALSE)
  }else{
  max_y = max(c(ZB_total))
  legend(x= (max_x + 2), y = max_y, legend = c("All Zooplankton", as.character(param$groups$species[zoo_groups])),
           col = c("Red",  coll), lty = c(1, rep(2, length(zoo_groups))), lwd = c(2, rep(1.5, length(zoo_groups))),bty = "n", cex = 0.8)
  par(xpd = FALSE)
  }
  }
}

## PLOT BIOMASS CONTRIBUTIONS
Bio_Cont_Plot = function(models, params){
modelss = models
param = params

w_min = which(round(log10(model$w), digits = 2) == -6.2)
w_max = which(round(log10(model$w), digits = 2) == -2)
N_sav = modelss$N
N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
N_ave = N_ave[,(w_min:w_max)]
zoo_groups =  which(is.na(param$groups$prop) == FALSE & rowSums(N_ave) != 0)
N_ave = N_ave[zoo_groups,]
Biom_ave = sweep(N_ave, 2, model$w[w_min:w_max], "*")
Biom_props = sweep(Biom_ave,2, colSums(Biom_ave), "/")

par(mar = c(4,4,4,8))
x = log10(model$w[w_min:w_max])
plot(x, rep(1,length(x)), type = "n", ylim = c(0.01,1), xlab = expression(paste("log"[10], "(Body Weight, g)")),
     ylab = "% Biomass Contribution")
coll = rainbow(length(zoo_groups))
for(i in 1:length(zoo_groups)){
  lines(x, Biom_props[i,], col = coll[i], lwd = 1.5, lty = 1)
}
title(main = paste("Phytoplankton Slope = " , round(enviro$b, digits = 2),
                   "\nTemperature = ", round(enviro$sst), "C"))
par(xpd = TRUE)
max_x = x[length(x)]
max_y = 1
legend(x= (max_x + 0.1), y = max_y, legend = c(as.character(param$groups$species[zoo_groups])),
       col = coll, lty = 1, lwd = 1.5,bty = "n", cex = 0.8)
par(xpd = FALSE)
}

## FISH DIET PLOTS
Fish_Diet = function(models, params){
  modelss = models
  param = params
  
  zoo_groups =  param$zoo_grps
  fish_groups = param$fish_grps
  num_zoo = param$num_zoo
  num_fish = param$num_fish
  
  ngrps = dim(param$groups)[1]
  w_min = which(log10(model$w) == -3)
  w_max = which(log10(model$w) == 6)
  assim_dynam = matrix(param$groups$alpha*param$nutrition, nrow = ngrps, ncol = ngrps, byrow = TRUE)/
                  matrix(param$nutrition, nrow = ngrps, ncol = ngrps)
  diet_sav = modelss$diet
  diet_ave = colMeans(diet_sav[(ceiling(0.5*dim(diet_sav)[1])):(dim(diet_sav)[1]),,,], dim = 1)
  diet_ave = diet_ave[,(w_min:w_max),]
  fish_ident = which(is.na(param$groups$prop) == TRUE)
  fish_diet = diet_ave[fish_ident,,]
  fish_diet_ingest = sweep(fish_diet, 2, c(1,assim_dynam[ngrps,]), "/")
  fish_diet_props = sweep(fish_diet_ingest, 1, rowSums(fish_diet_ingest), "/")
  fish_diet_props = fish_diet_props[,-1]
  zoo_groups =  which(apply(fish_diet_props, 2, max) > 0.03)
  zoo_groups = zoo_groups[which(is.na(param$groups$prop[zoo_groups]) == FALSE )]
  fish_diet_dom_zoo = fish_diet_props
  zoo_names = as.character(param$groups$species[1:num_zoo])  
  
  
  par(mar = c(4,4,4,2))
  x = log10(model$w[w_min:w_max])
  plot(x[-length(x)], rep(1, length(x[-length(x)])), type = "n", ylim = c(0, ceiling(max(fish_diet_dom_zoo[,1:9], na.rm = TRUE)*100)/100),
        xlab = expression(paste("log"[10], "(Body Weight, g)")), ylab = "% of Diet")
  
  coll = rainbow(9)
  for(i in 1:9){
    lines(x[-length(x)], fish_diet_dom_zoo[-length(x),i], lty = 1, lwd = 2, col = coll[i])
  }
  legend("topright", legend = zoo_names, lty = 1, col = coll,
         lwd = 2, bty = "n", cex = 0.8)
  title(main = paste("Phytoplankton Slope = " , round(enviro$b, digits = 2),
                     "\nTemperature = ", round(enviro$sst), "C"))
}

########## RESULTS TABLES
Summary_Results = function(modelss, param){
  N_sav = modelss$N
  N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
  zoo_groups =  param$zoo_grps
  fish_groups = param$fish_grps
  num_zoo = param$num_zoo
  num_fish = param$num_fish
  tot_zoo = colSums(N_ave[zoo_groups,])
  tot_fish = N_ave[fish_groups,]
  if(num_fish > 1){
    tot_fish = colSums(N_ave[fish_groups,])
  }
  
  ## BIOMASS OF DIFFERENT GROUPS
  zoo_bioms = rowSums(N_ave[zoo_groups,]*matrix(modelss$w, nrow = length(zoo_groups), 
                                                ncol = dim(N_ave)[2], byrow = TRUE))
  fish_bioms = rowSums(N_ave[fish_groups,]*matrix(modelss$w, nrow = length(fish_groups), 
                                    ncol = dim(N_ave)[2], byrow = TRUE))
 # fish_s = seq(1,which(round(log10(modelss$w), digits = 2) == 2),1) # (10^-3g - 10^2g)
#  fish_l = seq((length(fish_s)+1), length(modelss$w), 1) # (10^2g - 10^6g)
 # tot_small_fish = sum(N_ave[-zoo_groups, fish_s]*modelss$w[fish_s])
#  tot_large_fish =  sum(N_ave[-zoo_groups, fish_l]*modelss$w[fish_l])
 
  tot_bioms = matrix(c(zoo_bioms, fish_bioms), nrow = (num_zoo + num_fish), ncol = 1)
  rownames(tot_bioms) = as.character(param$groups$species)
  colnames(tot_bioms) = "Biomass g m^-3"
  tot_bioms = round(tot_bioms, digits = 4)
  
  ## Calculate information tables for zooplankton and fish communities
  FB_total = modelss$Biomass[(ceiling(0.5*dim(modelss$Biomass)[1])):(dim(modelss$Biomass)[1]), -c(zoo_groups)]
  Z_groups = modelss$Biomass[(ceiling(0.5*dim(modelss$Biomass)[1])):(dim(modelss$Biomass)[1]), zoo_groups]
  ZB_total = rowSums(Z_groups) 
  if(num_fish > 1){FB_total = rowSums(FB_total)}
  
  # Average Biomass
  ave_zoo = round(mean(ZB_total), digits = 2)
  ave_fish = round(mean(FB_total), digits = 2)
  ave_phyto = round(sum(model$nPP*model$w_phyto), digits = 2)
  
  # Average Slope
  fish_start = which(round(log10(model$w), digits = 2) == (param$groups$W0[dim(param$groups)[1]]))
  fish_finish = which(round(log10(model$w), digits = 2) == (param$groups$Wmat[dim(param$groups)[1]])) 
  max_phyto = round(log10(param$wMax_phyto), digits = 2)
  zoo_start =  which(round(log10(model$w), digits = 2) == -10.7) # 100um ESD
  zoo_finish = which(round(log10(model$w), digits = 2) == 0) # 1mm ESD
 zoo_slope2 = round(lm(log10(tot_zoo[zoo_start:zoo_finish])~log10(model$w[zoo_start:zoo_finish]))$
                       coefficients[2], digits = 2)
  #zoo_slope2 =  round((log10(tot_zoo[zoo_finish]) - log10(tot_zoo[zoo_start]))/
  #                      (log10(model$w[zoo_finish]) - log10(model$w[zoo_start])), digits = 2)
  fish_slope2 = round(lm(log10(tot_fish[fish_start:fish_finish])~log10(model$w[fish_start:fish_finish]))$
                        coefficients[2], digits = 2)
  #fish_slope2 = round((log10(N_ave[dim(N_ave)[1],fish_finish]) - log10(N_ave[dim(N_ave)[1],fish_start]))/
  #                      (log10(model$w[fish_finish]) - log10(model$w[fish_start])), digits = 2)
  phyto_slope = param$environ$b
  
  # Proportions and PPMR
  prop_zoo_biom = round(colMeans(Z_groups)/sum(colMeans(Z_groups)), digits = 2) # Biomass Proportions
  zoo_ave_m = round(sum(colMeans(Z_groups)/sum(colMeans(Z_groups))* # Average m-value
                          param$groups$m[zoo_groups]), digits = 2)
  
  NN = rowSums(modelss$N, dim = 2)
  
  prop_zoo_abund =  (colMeans(NN[(ceiling(0.5*dim(NN)[1])):(dim(NN)[1]), zoo_groups])/
                            sum(colMeans(NN[(ceiling(0.5*dim(NN)[1])):(dim(NN)[1]), zoo_groups])))# Abundance Proportions
  
  results = matrix(c(ave_phyto, ave_zoo, ave_fish, round(param$environ$b, digits = 2), zoo_slope2, fish_slope2,NA, zoo_ave_m, NA),
                    nrow = 3, ncol = 3)
  colnames(results) = c("Biomass", "Slope", "m")
  rownames(results) = c("Phytoplankton","Zooplankton", "Fish")
  
  zoo_results = matrix(c(prop_zoo_abund, prop_zoo_biom), nrow = num_zoo, ncol = 2)
  colnames(zoo_results) = c("% Abundance Total", "% Biomass Total")
  rownames(zoo_results) = as.character(param$groups$species[zoo_groups])
  
  return(list("General" = results, "Zoo_Groups" = zoo_results, "Total_Bioms" = tot_bioms))
}

Zoo_Abund_Pie = function(modelss, cut_point1, cut_point2){
  N_sav = modelss$N
  N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
  zoo_groups =  param$zoo_grps[-c(1,2)]
  num_zoo = param$num_zoo - 2
  weight_cut = which(modelss$w >= 10^cut_point1 & modelss$w <= 10^cut_point2)
  zoo_abunds = rowSums(N_ave[zoo_groups, weight_cut])
  
  par(mfrow = c(1,1))
  
  slices = zoo_abunds
  lbls = c("Larvaceans", "Omni Cop", "Carn Cop", "Euphausiids", "Chaetognaths", "Salps", "Jellyfish")
  pct = round(slices/sum(slices)*100)
  lbls = paste(lbls, pct)
  lbls = paste(lbls, "%", sep = "")
  pie(slices, labels = lbls, col = rainbow(length(lbls)), 
      main = "Zoo Abundances")
}

Zoo_Biom_Pie = function(modelss, cut_point1, cut_point2){
  N_sav = modelss$N
  N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
  zoo_groups =  param$zoo_grps[-c(1,2)]
  num_zoo = param$num_zoo - 2
  weight_cut = which(modelss$w >= 10^cut_point1 & modelss$w <= 10^cut_point2)
  zoo_bioms = sweep(N_ave[zoo_groups,], 2, modelss$w, "*")
  zoo_biomss = rowSums(zoo_bioms[, weight_cut])
  
  par(mfrow = c(1,1))
  slices = zoo_biomss
  #slices = c(0.21, 0.18, 0.14, 0.19, 0.19, 0.07, 0.02)
  lbls = c("Larvaceans", "Omni Cop", "Carn Cop", "Euphausiids", "Chaetognaths", "Salps", "Jellyfish")
  pct = round(slices/sum(slices)*100)
  lbls = paste(lbls, pct)
  lbls = paste(lbls, "%", sep = "")
  pie(slices, labels = lbls, col = rainbow(length(lbls)), 
      main = "", cex = 1.7)
} 

PPMR_plot = function(modelss, param){

  N_sav = modelss$N
  N_ave = colMeans(modelss$N[(ceiling(0.5*dim(N_sav)[1])):(dim(N_sav)[1]),,], dim = 1)
  large_zoo_bioms = sweep(N_ave[c(3:9),], 2, model$w, "*")
  large_zoo_m = param$groups$m
  large_zoo_m = large_zoo_m[3:9]
  
  D.z = 2*(3*(model$w)*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
  betas =  log10(t(sapply(large_zoo_m, function(x){(exp(0.02*log(D.z)^2 - x + 1.832))^3})))
  
  beta_props = large_zoo_bioms/sum(large_zoo_bioms)
  plot(density(betas, weights = beta_props), main = "", xlab = expression(paste("log"[10], " PPMR")),
       ylab = "Proportion of Zoo. Biomass", xlim = c(0, 10), 
       ylim = c(0,ceiling(max(density(betas, weights = beta_props)$y)*10)/10))
  abline(v = round(sum(beta_props*betas), digits = 2))
  title(main = paste("Phytoplankton Slope = " , round(enviro$b, digits = 2),
                     "\nTemperature = ", round(enviro$sst), "C",
                     "\n Ave. PPMR = ", round(sum(beta_props*betas), digits = 2)))
  
}



#groups = data.frame("species" = c("Flagellates", "Ciliates" ,   "Larvaceans",  "Omni Cop"  , 
#                                  "Carn Cop",    "Euphausiids", "Chaetog"  ,   "Salps" ,     
#                                  "Jellyfish",   "Fish_Small"  ,"Fish_Med"   , "Fish_Large" ),
#                    "prop" = c(1.000, 0.500, 0.040, 0.250, 0.350, 0.150, 0.250, 0.060, 0.005, NA,
#                               NA, NA),
#                    "W0" = c(-12.3,  -9.3,  -6.4 , -6.4  ,-6.4  ,-6.4  ,-5.0  ,-4.5,  -3.0,  -3.0,
#                             -3.0 , -3.0),
#                    "Wmax" = c( -6.3, -6.3, -2.3, -2.1, -1.2,  0.0,  0.0,  0.6,  2.3,  2.0,  4.0,  6.0),
#                    "Wmat" = c( -8.3, -8.3, -4.3, -4.1, -3.2,  -2, -2,  -1.4,  0.3,  0,  2.0,  4.0),
#                    "gamma" = as.integer(c(rep(7450, 9), rep(6400, 3))),
#                    "q" = c(rep(1, 9), rep(0.8, 3)),
#                    "m" = c(1.5, 0.04, -3, -0.48, 1.5, -0.5, 1, -2.68, 0.73, NA, NA, NA),
#                    "beta" = c(rep(NA, 9), rep(100, 3)),
#                    "sigma" = c(0.36, 0.47, 0.7, 0.57, 0.42, 0.57, 0.46, 0.7, 0.52, 1.3, 1.3, 1.3),
#                    "alpha" = c(rep(0.25, 9), rep(0.2, 3)),
#                    "type" = c("O", "O", "O", "H", "C", "H", "C", rep("O",5)),
#                    "carbon" = c(0.2, 0.2, 0.02, 0.12, 0.12, 0.12, 0.04, 0.02, 0.005, 0.1, 0.1, 0.1),
#                    "repro" = as.integer(rep(0,12)))
