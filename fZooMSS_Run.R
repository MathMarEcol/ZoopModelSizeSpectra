## Run model forward in time using previosly defined setup and model lists

## Last Updated Tuesday 17th March 2020
##

fZooMSS_Run <- function(model){

  # Pull out some useful parameters - just a shortcut
  param <- model$param
  dt <- model$param$dt
  dx <- model$param$dx
  ngrps <- model$param$ngrps
  ngrid <- model$param$ngrid
  w <- model$param$w
  W0 <- param$Groups$W0
  Wmax <- param$Groups$Wmax
  fish_grps <- model$param$fish_grps
  assim_eff <- model$assim_eff
  temp_eff <- model$temp_eff
  dynam_dietkernel <- model$dynam_dietkernel
  dynam_growthkernel <- model$dynam_growthkernel
  dynam_mortkernel <- model$dynam_mortkernel
  dynam_diffkernel <- model$dynam_diffkernel
  diff_phyto <- model$diff_phyto
  ingested_phyto <- model$ingested_phyto

  curr_min_size <- vector()
  curr_max_size <- vector()
  for (i in 1:ngrps){
    curr_min_size[i] <- which(round(log10(w), digits = 2) == W0[i])
    curr_max_size[i] <- which(round(log10(w), digits = 2) == Wmax[i])
  }

  idx_iter <- 2:ngrid
  idx <- 2:(ngrid-1)
  itimemax  <- param$tmax / dt  #max index of time array

  if(length(param$zoo_grps) > 1){ # If there's only one zoo group, then you do not need w0idx. All this stuff gives you info about all zoo groups except the smallest zoo group.
    w0idx <- which(W0 > min(W0) & is.na(param$Groups$Prop) == FALSE)
    w0mins <- rep(0, length(w0idx))
    props_z <- param$Groups$Prop[w0idx] # Zooplankton proportions

    for(i in 1:length(w0idx)){
      # Which size class is the smallest size class for each functional group
      w0mins[i] <- which(round(log10(w), digits = 2) == W0[w0idx[i]])
    }
  }

  # Matrices for MvF and MvF-D numeric solution
  A_iter <- matrix(0, nrow = ngrps, ncol = ngrid)
  C_iter <- matrix(0, nrow = ngrps, ncol = ngrid)
  S_iter <- matrix(0, nrow = ngrps, ncol = ngrid)

  A <- matrix(0, nrow = ngrps, ncol = ngrid)
  B <- matrix(0, nrow = ngrps, ncol = ngrid)
  C <- matrix(0, nrow = ngrps, ncol = ngrid)
  S <- matrix(0, nrow = ngrps, ncol = ngrid)

  # Temporary Matrices that get updated each time step some of these saved for output
  N <- matrix(model$N[1,,], nrow = ngrps, ncol = ngrid) # Abundances of functional groups, dim 1 = groups, dim 2 = size classes

  pb <- txtProgressBar(min = 0, max = itimemax, initial = 1, style = 3) # Initial progress bar

  # BIG TIME LOOP
  for (itime in 1:itimemax){

    setTxtProgressBar(pb, itime) # Update progress bar

    growth_multiplier <- colSums(N * assim_eff) # 1 x n_sizes
    predation_multiplier <- N * temp_eff # n_species x n_sizes
    diffusion_multiplier <- colSums(N * (assim_eff^2)) # 1 x n_sizes

    ### DO GROWTH
    dim(dynam_growthkernel) <- c(ngrps*ngrid, ngrid)
    cs <- .colSums(growth_multiplier * t(dynam_growthkernel), m = ngrid, n = ngrps*ngrid)
    dim(cs) <- c(ngrps, ngrid)
    gg <- ingested_phyto + cs

    # sw <- sweep(dynam_growthkernel, 3, growth_multiplier, '*') # n_species x n_sizes x n_sizes
    # ap <- aperm(sw, c(3,1,2)) # n_sizes x n_species x n_sizes
    # cs <- colSums(ap) # n_species x n_sizes
    # gg <- model$ingested_phyto + cs
    # rm(sw, ap, cs)

    ### DO MORTALITY
    # dynam_mortkernel <- model$dynam_mortkernel
    # cs <- colSums(predation_multiplier * t(dynam_mortkernel))
    # dim(cs) <- c(ngrps, ngrid)
    # gg <- ingested_phyto + cs

    sw2 <- sweep(dynam_mortkernel, c(2,3), predation_multiplier, '*') # n_sizes x n_species x n_sizes
    ap2 <- aperm(sw2, c(2,3,1))
    M2 <- .colSums(colSums(ap2),ngrid,ngrid) # 1 x n_sizes
    Z <- sweep(model$M_sb + model$fish_mort, 2, M2, '+') # Total dynamic spectrum mortality (n_species x n_sizes)
    rm(sw2, ap2)

    ### DO DIFFUSION
    dim(dynam_diffkernel) <- c(ngrps*ngrid, ngrid)
    cs <- .colSums(diffusion_multiplier * t(dynam_diffkernel), m = ngrid, n = ngrps*ngrid)
    dim(cs) <- c(ngrps, ngrid)
    diff <- diff_phyto + cs

    # sw3 <- sweep(dynam_diffkernel, 3, diffusion_multiplier, '*')
    # ap3 <- aperm(sw3, c(3,1,2))
    # cs3 <- colSums(ap3)
    # diff <- diff_phyto + cs3
    # rm(sw3, ap3, cs3)

    ### MvF WITH DIFFUSION ALGORITHM
    # Numerical implementation matrices (for MvF without diffusion)
    A_iter[,idx_iter] <- dt/dx * gg[,idx_iter-1] # Growth stuff
    C_iter[,idx_iter] <- 1 + dt * Z[,idx_iter] + dt/dx * gg[,idx_iter] # Mortality
    S_iter[,idx_iter] <- N[,idx_iter] # N at.....
    N_iter <- N # Current Abundance
    N_iter[1,1] <- N[1,1] # This forces R to make a copy of the variable. Otherwise N is linked to N_iter in the Rcpp code and they change together.

    # Numerical implementation matrices (for MvF WITH diffusion)
    A[,idx] <- dt/dx * (gg[,idx-1] + diff[,idx-1] * (log(10)/2+1/(2*dx))) # Growth stuff
    B[,idx] <- diff[,idx+1] * dt/(2*dx^2) # Diffusion term
    C[,idx] <- 1 + dt * Z[,idx] + dt/dx*(gg[,idx] + diff[,idx] * (log(10)/2+1/dx)) # Mortality
    S[,idx] <- N[,idx]

    # The original Base R code for the MvF equation
    N <- fZooMSS_MvF_BaseR(ngrps, curr_min_size, curr_max_size,
                            A_iter, C_iter, N_iter, S_iter,
                            A, B, C, N, S)

    # N <- fZooMSS_MvF_Rcpp(cngrps=ngrps, cN_iter=N_iter,
    #                        cA_iter=A_iter, cC_iter=C_iter, cS_iter=S_iter,
    #                        cN=N, cA=A, cB=B, cC=C, cS=S,
    #                        ccurr_min_size=curr_min_size, ccurr_max_size=curr_max_size)


    #### Keep smallest fish community size class as equal to equivalent zooplankton size class
    ### Keep smallest zooplankton size class abundnace for each group locked to others in size spectrum
    if(length(param$zoo_grps) > 1){ # If you only have one zoo group, it will be locked to phyto spectrum so you do not need to do this
      for(i in 1:length(w0idx)){
        w_min_curr <- w0mins[i]
        exclude_mins <- w0idx[which(w0mins == w_min_curr)]
        N[w0idx[i], w_min_curr] <- props_z[i] * sum(N[-exclude_mins, w_min_curr])
      }
    }

    fish_mins <- unlist(lapply(W0[fish_grps],
                               function(x){which(round(log10(w), digits = 2) == x)}))

    if(length(fish_grps) > 1 & length(param$zoo_grps) > 1){
      N[fish_grps,fish_mins] <- (1/length(fish_grps))*(colSums(N[-fish_grps,fish_mins]))
    }else{
      N[fish_grps, fish_mins] <- (1/length(fish_grps))*sum(N[-fish_grps, fish_mins])
    }


    # Save results:
    if((itime %% param$isave) == 0){
      isav <- itime/param$isave

      ## Phytoplankton diet
      pico_phyto_diet <- rowSums(model$diet_pico_phyto*N) # Pico-phytoplankton
      nano_phyto_diet <- rowSums(model$diet_nano_phyto*N) # Nano-phytoplankton
      micro_phyto_diet <- rowSums(model$diet_micro_phyto*N) # Micro-phytoplankton

      ## Functional group diet
      ### Create an ngrps*ngrid*ngrps*ngrid array of abundances, to save time without sweeps: dim1 = pred groups, dim 2 = pred sizes, dim 3 = prey groups, dim 4 = prey sizes
      N_array_temp <- aperm(replicate(ngrid, N), c(3,1,2))
      N_array <- aperm(replicate(ngrps, N_array_temp), c(4,1,2,3))
      dynam_diet <- rowSums(aperm(rowSums(sweep(dynam_dietkernel*N_array, c(1,2), N, "*"), dims = 3), c(1,3,2)), dims = 2)

      model$diet[isav,,1:3] <- cbind(pico_phyto_diet, nano_phyto_diet, micro_phyto_diet)
      model$diet[isav,,c(4:(dim(param$Groups)[1]+3))] <- dynam_diet
      model$N[isav,,] <- N # Save N by taxa and size
      model$Z[isav,,] <-  Z ## Save mortality
      model$gg[isav,,] <-  gg ## Save growth
    }
  } # End of time loop
  return(model)

}