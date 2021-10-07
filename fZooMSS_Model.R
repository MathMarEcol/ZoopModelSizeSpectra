fZooMSS_Model <- function(input_params, Groups, SaveTimeSteps = FALSE){

  source("fZooMSS_Params.R")
  source("fZooMSS_Setup.R")
  source("fZooMSS_MvF_BaseR.R")
  # sourceCpp("fZooMSS_MvF_Rcpp.cpp")
  source("fZooMSS_Run.R")
  source("fZooMSS_Xtras.R")

  input_params <- untibble(input_params)

  ################### RUN THE MODEL ###################
  param <- fZooMSS_Params(Groups, input_params) # Set up parameter list
  model <- fZooMSS_Setup(param) # Set up model equation stuff
  model_output <- fZooMSS_Run(model) # Run the model

  ################### AVERAGE THE LAST 50 % OF THE MODEL RUN ###################

  ave_abundances <- fZooMSS_AveOutput(model_output$N)
  ave_diets <- fZooMSS_AveOutput(model_output$diet)
  ave_growth <- fZooMSS_AveOutput(model_output$gg)
  ave_mort <- fZooMSS_AveOutput(model_output$Z)

  if (SaveTimeSteps == TRUE){
    model_output$Abundance <- rowSums(model$N, dims = 2) ## Save Total Abundance
    model_output$Biomass <- colSums(aperm(sweep(model$N, 3, model$param$w, "*"), c(3,1,2)))

  results <- list("abundances" = ave_abundances, # Save mean abundance
                 "diets" = ave_diets,  # Save mean diets
                 "growth" = ave_growth,  # Save mean growth
                 "mortality" = ave_mort, # Save mean predation
                 "model" = model_output) # Save whole model
  }
  if (SaveTimeSteps == FALSE){
    reduce_output <- list(param = model_output$param) # Create a new list so we can have identical structure of model$param
    results <- list("abundances" = ave_abundances, # Save mean abundance
                   "diets" = ave_diets,  # Save mean diets
                   "growth" = ave_growth,  # Save mean growth
                   "mortality" = ave_mort, # Save mean predation
                   "model" = reduce_output) # Save parameters only
  }

  return(results)
}
