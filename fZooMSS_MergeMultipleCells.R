# This script takes a run name and loads all
# the data from the individual cells and
# merges them into a single list
#
#
# Written by Jason Everett (UQ/CSIRO/UNSW)
# Written: Sunday 2nd February 2020
# Last Updated: Friday 9th October 2020

fZooMSS_MergeMultipleCells <- function(){
  run <- basename(getwd())

  library(tidyverse)
  files <- sort(list.files("RawOutput", full.names = TRUE))

  res <- list()
  diets <- list()
  growth <- list()

  for (f in 1:length(files)){
    out <- read_rds(files[f])

    res[[f]] <- out$abundances
    growth[[f]] <- out$growth
    diets[[f]] <- out$diets

    if (f == 1){
      mdl <- out$model
    }

    rm(out)
  }

  saveRDS(mdl, paste0("Output/model_",run,".RDS"))
  saveRDS(res, paste0("Output/res_",run,".RDS"))
  saveRDS(growth, paste0("Output/growth_",run,".RDS"))
  saveRDS(diets, paste0("Output/diets_",run,".RDS"))

  save(list = c("res", "growth", "diets", "mdl"), file = paste0(paste0("Output/full_",run,".RData")))

  print("Consider running `tar -cf RawOutput.tar RawOutput` to reduce the number of files to be synced")

}
