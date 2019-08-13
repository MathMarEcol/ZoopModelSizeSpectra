setwd("~/Desktop/Multi Core Code")
Groups <- read.csv("Test Groups.csv")
source("Multi-Zoo Model.R")

enviro_vector <- enviro_data[1,]
param <- params(Groups, enviro_vector, tmax = 1, f_mort = 0) # Set up parameter list
model <- Setup(param) # Set up model equation stuff 

slurm_fix <- function(slr_job, outtype){
  if (!(class(slr_job) == "slurm_job")) {
    stop("slr_job must be a slurm_job")
  }
  outtypes <- c("table", "raw")
  if (!(outtype %in% outtypes)) {
    stop(paste("outtype should be one of:", paste(outtypes, 
                                                  collapse = ", ")))
  }
  res_files <- paste0("results_", 0:(slr_job$nodes + 00 - 1), ".RDS")
  tmpdir <- paste0("_rslurm_", slr_job$jobname)
  missing_files <- setdiff(res_files, dir(path = tmpdir))
  if (length(missing_files) > 0) {
    missing_list <- paste(missing_files, collapse = ", ")
    warning(paste("The following files are missing:", missing_list))
  }
  res_files <- file.path(tmpdir, setdiff(res_files, missing_files))
  if (length(res_files) == 0) 
    return(NA)
  slurm_out <- lapply(res_files, readRDS)
  slurm_out <- do.call(c, slurm_out)
  if (outtype == "table") {
    slurm_out <- as.data.frame(do.call(rbind, slurm_out))
  }
  slurm_out
}

######## SAVE ITERATION
save_iteration <- function(filename, pics){
  
  new_file <- paste("./Iterations/", filename, "/", sep = "")
  dir.create(new_file)
  
  ### SAVE f (model code), ADD_OBJECTS (Test Groups), PARAMS (enviro_data)
  file.copy('./_rslurm_multi_zoo_model/f.RDS', new_file)
  file.copy('./_rslurm_multi_zoo_model/add_objects.RData', new_file)
  file.copy('./_rslurm_multi_zoo_model/params.RDS', new_file)
  
  ### SAVE GROWTH, DIETS AND ABUNDANCES
  save(growth, file = paste(new_file, "growth.RData", sep = "")) # Growth
  save(res, file = paste(new_file, "res.RData", sep = "")) # Abundances
  save(diets, file = paste(new_file, "diets.RData", sep = "")) # Diets
  
  ### SAVE PLOTS OF BIOM PROPS AGAINST CHLO
  if(pics == TRUE){
    filename = paste(new_file,"BiomProps_Chlo", ".pdf" , sep = "")
    pdf(filename, width = 5.8, height = 8.3)
    biom_props(res, groups = Groups, -10.7, 3, en = "chlo")
    dev.off()
  }
}



#### PLOT AVERAGE SPECTRUM
Spectrum_Plot <- function(N, groups, i){
  N_ave = N[[i]]
  enviro <- enviro_data[i,]
  fish_groups = which(is.na(groups$prop) == TRUE) # Which rows are fish
  zoo_groups = which(is.na(groups$prop) == FALSE) # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  num_fish = length(fish_groups)
  nPP = 10^(enviro$a)*(model$w_phyto^(enviro$b))
  tot_zoo = colSums(N_ave[zoo_groups,])
  
  if(num_fish > 0){
    if(num_fish > 1){
      tot_fish = colSums(N_ave[fish_groups,])
    }else{tot_fish = N_ave[fish_groups,]}
  }
  
  y_d_ref = which(round(log10(model$w), digits = 2) == 5)
    y_down = floor(log10(N_ave[dim(N_ave)[1],][y_d_ref]))
  y_up = ceiling(log10(model$nPP[5]))
  
  ## PLOT PHYTO-ZOO-FISH TOTALS
  par(mfrow = c(1,1), mar = c(4,4,4,2))
  plot(log10(model$w_phyto), log10(nPP), type= "l", col = "green",
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
  coll = rainbow(num_zoo)
  
  for(i in 1:num_zoo){
    colll = coll[i]
    lines(log10(model$w), log10(N_ave[zoo_groups[i],]), lty = 2, col = colll, lwd = 1.5)
  }
  legend(x = 0, y = y_up + 1.5, legend = as.character(groups$species[zoo_groups]),
         col = coll, lty = 2, lwd = 1.5,bty = "n", cex = 0.8)
  
  ## PLOT FISH GROUPS
  if(num_fish > 0){
    for(i in 1:num_fish){
      lines(log10(model$w), log10(N_ave[fish_groups[i],]), lty = 2, col = "blue", lwd = 1.5)
    }
  }
}

####### SLOPES
slope_plots <- function(N, groups, start, finish){
  ########## RESULTS TABLES
  slopes <- matrix(0, dim(enviro_data)[1],3)
  fish_groups = which(is.na(groups$prop) == TRUE) # Which rows are fish
  zoo_groups = which(is.na(groups$prop) == FALSE) # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  num_fish = length(fish_groups)
  
  for(i in 1:dim(enviro_data)[1]){
  N_ave = N[[i]]
  enviro <- enviro_data[i,]
  tot_zoo = colSums(N_ave[zoo_groups,])
  tot_fish = N_ave[fish_groups,]
  
    if(num_fish > 1){
      tot_fish = colSums(N_ave[fish_groups,])
    }
    tot_anim = tot_zoo + tot_fish
    # Average Slope
    fish_start = which(round(log10(model$w), digits = 2) == 0)
    fish_finish = which(round(log10(model$w), digits = 2) == (param$groups$Wmat[dim(groups)[1]])) 
    max_phyto = round(log10(param$wMax_phyto), digits = 2)
    zoo_start =  which(round(log10(model$w), digits = 2) == start)
    zoo_finish = which(round(log10(model$w), digits = 2) == finish) 
    zoo_slope2 = lm(log10(tot_anim[zoo_start:zoo_finish])~log10(model$w[zoo_start:zoo_finish]))$coefficients[2]
    #zoo_slope2 =  round((log10(tot_zoo[zoo_finish]) - log10(tot_zoo[zoo_start]))/
    #                      (log10(model$w[zoo_finish]) - log10(model$w[zoo_start])), digits = 2)
    fish_slope2 = round(lm(log10(tot_fish[fish_start:fish_finish])~log10(model$w[fish_start:fish_finish]))$
                          coefficients[2], digits = 2)
    #fish_slope2 = round((log10(N_ave[dim(N_ave)[1],fish_finish]) - log10(N_ave[dim(N_ave)[1],fish_start]))/
    #                      (log10(model$w[fish_finish]) - log10(model$w[fish_start])), digits = 2)
    phyto_slope = enviro$b
    
    slopes[i,1] <- phyto_slope
    slopes[i,2] <- zoo_slope2
    slopes[i,3] <- fish_slope2
  }
  
  slopes[slopes[,2] > -0.85, 2] <- NA
  
  par(mfrow = c(3,1))
  plot(log10(enviro_data$chlo), slopes[,1], ylab = "Slope", xlab = "log10(Chlo)", main = "Phyto Slope")
  plot(log10(enviro_data$chlo), slopes[,2], ylab = "Slope", xlab = "log10(Chlo)", main = "Zoo Slope")
  plot(log10(enviro_data$chlo), slopes[,2] - slopes[,1], ylab = "Zoo Slope - Phyto Slope", xlab = "log10(Chlo)", main = "Slope Diff")
  abline(0,0, col = "red", lwd = 2)
  par(mfrow = c(1,1))
  
}

####### PROPORTIONS OF ABUNDANCE
abund_props <- function(N, groups, cut_point1, cut_point2, en){
  zoo_groups = which(is.na(groups$prop) == FALSE)[-c(1,2)] # Which rows are zooplankton
  num_zoo = length(zoo_groups)

  abund_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)

  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    tot_zoo = colSums(N_ave[zoo_groups,])

   weight_cut = which(round(log10(model$w),2) >= cut_point1 & round(log10(model$w),2) <= cut_point2)
   zoo_abunds = rowSums(N_ave[zoo_groups, weight_cut])
   zoo_props = zoo_abunds/sum(zoo_abunds)
   abund_mat[i,] <- zoo_props 
  }
  

  if(en == "chlo"){
    par(mfrow = c(4,2))
    for( i in 1:num_zoo){
      plot(log10(enviro_data[,"chlo"]), abund_mat[,i], main = groups$species[i+2], xlab = "log10(Chlo)", ylab = "Abund Prop")
    }
  }

  
  if(en == "sst"){
    par(mfrow = c(4,2))
    for( i in 1:num_zoo){
      plot((enviro_data$sst), abund_mat[,i], main = groups$species[i+2], xlab = "SST", ylab = "Abund Prop")
    }
  }
  par(mfrow = c(1,1))
  
  
}


####### ACTUAL ABUNDANCE
abunds_act <- function(N, groups, cut_point1, cut_point2, en){
  zoo_groups = which(is.na(groups$prop) == FALSE)[-c(1,2)] # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  
  abund_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)
  
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    tot_zoo = colSums(N_ave[zoo_groups,])
    
    weight_cut = which(round(log10(model$w),2) >= cut_point1 & round(log10(model$w),2) <= cut_point2)
    zoo_abunds = rowSums(N_ave[zoo_groups, weight_cut])
    zoo_props = zoo_abunds
    abund_mat[i,] <- zoo_props 
  }
  
  
  if(en == "chlo"){
    par(mfrow = c(4,2))
    for( i in 1:num_zoo){
      plot(log10(enviro_data[,"chlo"]), log10(abund_mat[,i]), main = groups$species[i+2], xlab = "log10(Chlo)", ylab = expression(paste("Abundance (# m"^-3, ")")))
    }
  }
  
  
  if(en == "sst"){
    par(mfrow = c(4,2))
    for( i in 1:num_zoo){
      plot((enviro_data$sst), log10(abund_mat[,i]), main = groups$species[i+2], xlab = "SST", ylab = expression(paste("Abundance (# m"^-3, ")")))
    }
  }
  par(mfrow = c(1,1))
  
  
}

####### PROPORTIONS OF BIOMASS
biom_props <- function(N, groups, cut_point1, cut_point2, en){
  zoo_groups = which(is.na(groups$prop) == FALSE)[-c(1,2)] # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  
  abund_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)
  
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    weight_cut = which(round(log10(model$w),2) >= cut_point1 & round(log10(model$w),2) <= cut_point2)
    zoo_abunds = rowSums(B_ave[zoo_groups, weight_cut])
    zoo_props = zoo_abunds/sum(zoo_abunds)
    abund_mat[i,] <- zoo_props 
  }
  
  if(en == "chlo"){
  par(mfrow = c(4,2))
  for( i in 1:num_zoo){
    plot(log10(enviro_data$chlo), abund_mat[,i], main = groups$species[i+2], xlab = "log10(Chlo)", ylab = "Biom Prop")
  }
  }
  
  if(en == "sst"){
    par(mfrow = c(4,2))
    for( i in 1:num_zoo){
      plot((enviro_data$sst), abund_mat[,i], main = groups$species[i+2], xlab = "SST", ylab = "Biom Prop")
    }
  }
  par(mfrow = c(1,1))
  
  
}

### ACTUAL BIOMASS
bioms_act <- function(N, groups, cut_point1, cut_point2, en){
  zoo_groups = which(is.na(groups$prop) == FALSE)[-c(1,2)] # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  
  abund_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)
  
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    weight_cut = which(round(log10(model$w),2) >= cut_point1 & round(log10(model$w),2) <= cut_point2)
    zoo_abunds = rowSums(B_ave[zoo_groups, weight_cut])
    zoo_props = zoo_abunds
    abund_mat[i,] <- zoo_props 
  }
  
  if(en == "chlo"){
    par(mfrow = c(4,2))
    for( i in 1:num_zoo){
      plot(log10(enviro_data$chlo), abund_mat[,i], main = groups$species[i+2], xlab = "log10(Chlo)", ylab = expression(paste("Biomass (g m"^-3, ")")))
    }
  }
  
  if(en == "sst"){
    par(mfrow = c(4,2))
    for( i in 1:num_zoo){
      plot((enviro_data$sst), abund_mat[,i], main = groups$species[i+2], xlab = "SST", ylab = expression(paste("Biomass (# m"^-3, ")")))
    }
  }
  par(mfrow = c(1,1))
  
  
}

### ABUND DATAFRAME
abund_mat <- function(N, groups, cut_point1, cut_point2, en){
  zoo_groups = which(is.na(groups$prop) == FALSE)[-c(1,2)] # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  
  abund_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)
  
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    tot_zoo = colSums(N_ave[zoo_groups,])
    
    weight_cut = which(round(log10(model$w),2) >= cut_point1 & round(log10(model$w),2) <= cut_point2)
    zoo_abunds = rowSums(N_ave[zoo_groups, weight_cut])
    zoo_props = zoo_abunds
    abund_mat[i,] <- zoo_props 
  }

  colnames(abund_mat) = as.character(groups$species)[3:9]
  abund_mat = as.data.frame(abund_mat)
  abund_mat
}


### BIOM DATAFRAME
biom_mat <- function(N, groups, cut_point1, cut_point2, en){
  zoo_groups = which(is.na(groups$prop) == FALSE)[-c(1,2)] # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  
  biom_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)
  
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    weight_cut = which(round(log10(model$w),2) >= cut_point1 & round(log10(model$w),2) <= cut_point2)
    zoo_bioms = rowSums(B_ave[zoo_groups, weight_cut])
    zoo_props = zoo_bioms
    biom_mat[i,] <- zoo_props 
  }
  
  colnames(biom_mat) = as.character(groups$species)[3:9]
  biom_mat = as.data.frame(biom_mat)
  biom_mat
}

### BIOM DATAFRAME
biom_mat_small <- function(N, groups, cut_point1, cut_point2, en){
  zoo_groups = which(is.na(groups$prop) == FALSE)[c(1,2)] # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  
  biom_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)
  
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    weight_cut = which(round(log10(model$w),2) >= cut_point1 & round(log10(model$w),2) <= cut_point2)
    zoo_bioms = rowSums(B_ave[zoo_groups, weight_cut])
    zoo_props = zoo_bioms
    biom_mat[i,] <- zoo_props 
  }
  
  colnames(biom_mat) = as.character(groups$species)[1:2]
  biom_mat = as.data.frame(biom_mat)
  biom_mat
}

slope_mat <- function(N, groups, start, finish){
  ########## RESULTS TABLES
  slopes <- matrix(0, dim(enviro_data)[1],3)
  fish_groups = which(is.na(groups$prop) == TRUE) # Which rows are fish
  zoo_groups = which(is.na(groups$prop) == FALSE) # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  num_fish = length(fish_groups)
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    enviro <- enviro_data[i,]
    tot_zoo = colSums(N_ave[zoo_groups,])
    tot_fish = N_ave[fish_groups,]
    
    if(num_fish > 1){
      tot_fish = colSums(N_ave[fish_groups,])
    }
    
    tot_anim = tot_zoo + tot_fish
    # Average Slope
    fish_start = which(round(log10(model$w), digits = 2) == (param$groups$W0[dim(groups)[1]]))
    fish_finish = which(round(log10(model$w), digits = 2) == (param$groups$Wmat[dim(groups)[1]])) 
    max_phyto = round(log10(param$wMax_phyto), digits = 2)
    zoo_start =  which(round(log10(model$w), digits = 2) == start)
    zoo_finish = which(round(log10(model$w), digits = 2) == finish) 
    zoo_slope2 = lm(log10(tot_anim[zoo_start:zoo_finish])~log10(model$w[zoo_start:zoo_finish]))$
      coefficients[2]
    #zoo_slope2 =  round((log10(tot_zoo[zoo_finish]) - log10(tot_zoo[zoo_start]))/
    #                      (log10(model$w[zoo_finish]) - log10(model$w[zoo_start])), digits = 2)
    fish_slope2 = round(lm(log10(tot_fish[fish_start:fish_finish])~log10(model$w[fish_start:fish_finish]))$
                          coefficients[2], digits = 2)
    #fish_slope2 = round((log10(N_ave[dim(N_ave)[1],fish_finish]) - log10(N_ave[dim(N_ave)[1],fish_start]))/
    #                      (log10(model$w[fish_finish]) - log10(model$w[fish_start])), digits = 2)
    phyto_slope = enviro$b
    
    slopes[i,1] <- phyto_slope
    slopes[i,2] <- zoo_slope2
    slopes[i,3] <- fish_slope2
  }
  
  colnames(slopes) = c("Phyto Slope", "Zoo Slope", "Fish Slope")
  slopes = as.data.frame(slopes)
  return(slopes)
  
}


### RATIOS
ratio_plot <- function(N, groups, enviro, cut_point1, cut_point2){
  enviro_data = enviro
  
  phyto_carb = rep(0, dim(enviro_data)[1])
  
  for(i in 1:dim(enviro_data)[1]){
  w0_phyto = -14.5		# minimum phytoplankton size class (1um)
  wMax_phyto = min(enviro_data$phyto_max[i])
  wMax_phyto = min(-8, wMax_phyto)
  w_phyto = 10^(seq(-14.5, wMax_phyto,0.1))
  a = enviro_data$a[i]
  b = enviro_data$b[i]
  phyto_carb[i] = sum(0.15*(10^a)*(w_phyto^(b+1)))
 # phyto_carb[i] = 0.15*((10^a/(b+1))*((10^wMax_phyto)^(b+1) - (10^w0_phyto)^(b+1)))
  }

  phyto_ww <- phyto_carb/0.15

  phyto_www <- 60*enviro_data$chlo/(0.1*1000)
    
  zoo_groups = which(is.na(groups$prop) == FALSE)
#  zoo_groups = c(1,2,4,5,6,7,8,9)# Which rows are zooplankton
  num_zoo = length(zoo_groups)
  
  biom_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol =12)
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    weight_cutz = c(round(log10(model$w),2) >=-10.7 & round(log10(model$w),2) <=3)
    weight_cutf = c(round(log10(model$w),2) >=-3 & round(log10(model$w),2) <= -3.1)
    weight_cut_mat = matrix(c(rep(weight_cutz, 9), rep(weight_cutf,3)), 
                            nrow = 12, ncol = length(model$w),byrow = TRUE)
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    zoo_bioms = rowSums(B_ave*weight_cut_mat)
    zoo_props = zoo_bioms
    biom_mat[i,] <- zoo_props 
  }
  
  #colnames(biom_mat) = as.character(groups$species)[zoo_groups]
  biom_mat = as.data.frame(biom_mat)
  zoo_ww = rowSums(biom_mat)
  zoo_carb = rowSums(sweep(biom_mat, 2, groups$carbon, "*"))
  
  fish_groups = which(is.na(groups$prop) == TRUE) # Which rows are fish
  num_fish = length(fish_groups)
  
  biom_matf <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_fish)
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    weight_cut = which(round(log10(model$w),2) >=-3 & round(log10(model$w),2) <= 6)
    fish_bioms = rowSums(B_ave[fish_groups,weight_cut])
    fish_props = fish_bioms
    biom_matf[i,] <- fish_props 
  }
  
  colnames(biom_matf) = as.character(groups$species)[10:12]
  biom_matf = as.data.frame(biom_matf)
  fish_ww = rowSums(biom_matf)
  fish_carb = rowSums(sweep(biom_matf, 2, groups$carbon[10:12], "*"))
  
  #phyto_ww = phyto_ww[-291]
  #zoo_ww = zoo_ww[-291]
  #fish_ww = fish_ww[-291]
  #enviro_data = enviro_data[-291,]
  par(mfrow = c(3,2))
    plot(log10(enviro_data$chlo), (phyto_ww), xlab = expression(paste("log"[10], "(Chlo)")),
         ylab = "", main = "Phyto Biomass")
    plot(log10(enviro_data$chlo), (zoo_ww/phyto_ww), xlab = expression(paste("log"[10], "(Chlo)")),
         ylab = "", main = "Zoo : Phyto")
    plot(log10(enviro_data$chlo), (zoo_ww), xlab = expression(paste("log"[10], "(Chlo)")),
         ylab = "", main = "Zoo Biomass")
 #   cols = rainbow(4)
#    biom_mat[,10] <- rowSums(biom_mat[,c(10:12)])
 #   choose <- c(3,4,6,8)
  #  for(i in 1:4){
  #     points(log10(enviro_data$chlo), biom_mat[,choose[i]], col= cols[i])
  #  }
    plot(log10(enviro_data$chlo), (fish_ww/zoo_ww), xlab = expression(paste("log"[10], "(Chlo)")),
         ylab = "", main = "Fish : Zoo")
    plot(log10(enviro_data$chlo), (fish_ww), xlab = expression(paste("log"[10], "(Chlo)")),
         ylab = "", main = "Fish Biomass")
    plot(log10(enviro_data$chlo), (fish_ww/phyto_ww), xlab = expression(paste("log"[10], "(Chlo)")),
         ylab = "", main = "Fish : Phyto")
    par(mfrow = c(1,1))
    
  }


### ZOO PPMR and CARBON VS CHLO
chlo_ppmr_carb <- function(N, groups, enviro, cut_point1, cut_point2){
  chlo <- enviro$chlo

  zoo_groups = which(is.na(groups$prop) == FALSE) # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  
  ave_zoo_ppmr <- rep(0, dim(enviro_data)[1])
  ave_zoo_carbon <- rep(0, dim(enviro_data)[1])
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    
    ### AVE ZOO PPMR
    weight_cut = which(round(log10(model$w),2) >= cut_point1 & round(log10(model$w),2) <= cut_point2)
    large_zoo_bioms = sweep(N_ave[c(3:9),], 2, model$w, "*")
    large_zoo_m = groups$m[3:9]
    
    D.z <- 2*(3*(model$w)*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
    betas =  log10(t(sapply(large_zoo_m, function(x){(exp(0.02*log(D.z)^2 - x + 1.832))^3})))
    
    beta_props = large_zoo_bioms[,weight_cut]/sum(large_zoo_bioms[,weight_cut])
    ave_zoo_ppmr[i] = (sum(beta_props*betas[,weight_cut]))
    
    ### AVE ZOO CARBON
    large_zoo_bioms = sweep(N_ave[c(3:9),], 2, model$w, "*")
    large_zoo_carb = sweep(large_zoo_bioms, 1, groups$carbon[3:9], "*")
    ave_zoo_carbon[i] = sum(large_zoo_carb[,weight_cut])/sum(large_zoo_bioms[,weight_cut])
  }
  
  par(mfrow = c(3,1))
  plot(log10(enviro_data$chlo), ave_zoo_ppmr, xlab = expression(paste("log"[10], "(Chlo)")),
       ylab = "", main = "Ave Zoo PPMR")
  plot(log10(enviro_data$chlo), ave_zoo_carbon, xlab = expression(paste("log"[10], "(Chlo)")),
       ylab = "", main = "Ave Zoo Proportion Carbon")
  plot(log10(enviro_data$chlo), ave_zoo_ppmr*ave_zoo_carbon, xlab = expression(paste("log"[10], "(Chlo)")),
       ylab = "", main = "Ave Zoo PPMR * Carbon")
  par(mfrow = c(1,1))
}

### PHYTO TO SMALL HETEROS
hp_plot <- function(N, enviro, en){
  
  hp_ratio = rep(0, dim(enviro)[1])
  for(i in 1:dim(enviro)[1]){
  environ = enviro[i,]
  
  ## How much phytoplankton is there
  phyto_biom <- sum(10^(environ$a)*(model$w_phyto^(environ$b+1)))
  
  ## How much hetero flag and cil is there
  N_curr = N[[i]]
  hetero_biom = sum(sweep(N_curr[c(1,2),],2, model$w, "*"))

  hp_ratio[i] = hetero_biom/phyto_biom
  }
  
  if(en == "chlo"){
  par(mfrow = c(1,1))
  plot(log10(enviro$chlo), hp_ratio, xlab = "log10(Chlo)", ylab = "Ratio", main = "Hetero F and C \n to Phyto Biomass")
  }
  
  if(en == "sst"){
    par(mfrow = c(1,1))
    plot((enviro$sst), hp_ratio, xlab = "SST", ylab = "Ratio", main = "Hetero F and C \n to Phyto Biomass")
  }
  
}

omni_plots <- function(diets, enviro){
  diet_store = matrix(0, nrow = dim(enviro)[1], ncol = 4)
  colnames(diet_store) = c("Ocops", "Larvs", "Krill", "Salps")
  zoo_curr = c(4,3,8,6)
  
  for(i in 1:4){
   for(j in 1:dim(enviro)[1]){
    group = zoo_curr[i]
    curr_diet = diets[[j]][group,-c(11:13)]
    diet_frac =  curr_diet/sum(curr_diet)
    diet_store[j,i] = sum(diet_frac[1:3])/sum(diet_frac)
   }
  }
  
  par(mfrow = c(2,2), oma = c(0, 0, 2, 0))
  for(i in 1:4){
    title("Phyto in Diet", outer = TRUE)
    plot(log10(enviro$chlo), diet_store[,i], xlab = "log10(Chlo)", ylab = "Proportion", ylim = c(0.5, 1),main = colnames(diet_store)[i])
  }
  par(mfrow = c(1,1))
  
}
  
### DIETS
diet_plot <- function(diets, group, enviro, en){
  diet_store = matrix(0, nrow = dim(enviro)[1], ncol = 15)
  colnames(diet_store) = c("SPhyto", "MPhyto", "LPhyto","Flag", "Cil", "Larvs", "Ocops", "Ccops", "Krill", "Chaets", "Salps", "Jelly", "Fish (<100g)", "Fish (100g-10kg)", "Fish (>10kg)")
  
  for(i in 1:dim(enviro)[1]){
  curr_diet = diets[[i]]
  curr_diet = curr_diet[group,]
  diet_store[i,] =  curr_diet/sum(curr_diet)
  }
  
  num_zoo = (which(colSums(diet_store)/sum(colSums(diet_store)) > 0.001))
  num_plot_rows = floor(length(num_zoo)/2 + length(num_zoo) %% 2)
  
  if(en == "chlo"){
    par(mfrow = c(num_plot_rows,2), oma = c(0,0,2,0))
    if(length(num_zoo) == 1){
      par(mfrow = c(1,1), oma = c(0,0,2,0))
    }
    
    for(i in 1:length(num_zoo)){
      curr_name = names(num_zoo)[i]
      if(group == 13){
      title(paste("<1 kg Fish", " Diet", sep = ""), outer = TRUE)
      }
      if(group != 13){
        title(paste(Groups$species[group], " Diet", sep = ""), outer = TRUE)
      }
      plot(log10(enviro_data$chlo), diet_store[,num_zoo[i]], main = curr_name, xlab = "log10(Chlo)", ylab = "Proportion")
    }
  }
  
  if(en == "sst"){
    par(mfrow = c(num_plot_rows,2), oma = c(0,0,2,0))
    if(length(num_zoo) == 1){
      par(mfrow = c(1,1), oma = c(0,0,2,0))
    }
    for(i in 1:length(num_zoo)){
      title(paste(Groups$species[group], " Diet", sep = ""), outer = TRUE)
      plot((enviro_data$sst), diet_store[,num_zoo[i]], main = c(names(num_zoo)[i]), xlab = "log10(Chlo)", ylab = "Proportion")
    }
  }
  par(mfrow = c(1,1))
}

pp_eat =matrix(0, nrow = 2, ncol = 200)
for(i in 1:dim(enviro_data)[1]){
  chlo = enviro_data$chlo[i]
  c_chl <- ((chlo^0.89)*(10^1.79))/chlo # chlo:carbon ratio (0.02 Chl:C)
  phyto_carb <- c_chl*chlo/1000 # (convert to grams carbon)
  phyto_prod <- phyto_carb/0.15 #(grams of wet weight)
  
  tot_pp_eat = sum(diets[[i]][,1])*50/365
  pp_eat[1,i] = phyto_prod
  pp_eat[2,i] = tot_pp_eat/phyto_prod
}

### RATIOS DATA FRAME
ratio_frame <- function(N, groups, enviro, cut_point1, cut_point2){
  enviro_data = enviro
  chlo <- enviro_data$chlo
  c_chl <- ((chlo^0.89)*(10^1.79))/chlo # chlo:carbon ratio (0.02 Chl:C)
  phyto_carb <- c_chl*chlo/1000 # (convert to grams carbon)
  phyto_ww <- phyto_carb*0.1
  
  
  zoo_groups = which(is.na(groups$prop) == FALSE)[-c(1:2)] # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  
  biom_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_zoo)
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    weight_cut = which(round(log10(model$w),2) >= cut_point1 & round(log10(model$w),2) <= cut_point2)
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    zoo_bioms = rowSums(B_ave[zoo_groups, weight_cut])
    zoo_props = zoo_bioms
    biom_mat[i,] <- zoo_props 
  }
  
  colnames(biom_mat) = as.character(groups$species)[3:9]
  biom_mat = as.data.frame(biom_mat)
  zoo_ww = rowSums(biom_mat)
  zoo_carb = rowSums(sweep(biom_mat, 2, groups$carbon[3:9], "*"))
  
  fish_groups = which(is.na(groups$prop) == TRUE) # Which rows are fish
  num_fish = length(fish_groups)
  
  biom_matf <- matrix(0, nrow = dim(enviro_data)[1], ncol = num_fish)
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    fish_bioms = rowSums(B_ave[fish_groups,])
    fish_props = fish_bioms
    biom_matf[i,] <- fish_props 
  }
  
  colnames(biom_matf) = as.character(groups$species)[10:12]
  biom_matf = as.data.frame(biom_matf)
  fish_ww = rowSums(biom_matf)
  fish_carb = rowSums(sweep(biom_matf, 2, groups$carbon[10:12], "*"))
  
  ratio_frame = data.frame("ZP" = zoo_carb/phyto_carb, "FZ" = fish_carb/zoo_carb, "FP" = fish_carb/phyto_carb)
  ratio_frame
  
}

slope_frame <- function(N, groups, start, finish){
  ########## RESULTS TABLES
  slopes <- matrix(0, dim(enviro_data)[1],3)
  fish_groups = which(is.na(groups$prop) == TRUE) # Which rows are fish
  zoo_groups = which(is.na(groups$prop) == FALSE) # Which rows are zooplankton
  num_zoo = length(zoo_groups)
  num_fish = length(fish_groups)
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    enviro <- enviro_data[i,]
    tot_zoo = colSums(N_ave[zoo_groups,])
    tot_fish = N_ave[fish_groups,]
    
    if(num_fish > 1){
      tot_fish = colSums(N_ave[fish_groups,])
    }
    tot_anim = tot_zoo + tot_fish
    # Average Slope
    fish_start = which(round(log10(model$w), digits = 2) == (param$groups$W0[dim(groups)[1]]))
    fish_finish = which(round(log10(model$w), digits = 2) == (param$groups$Wmat[dim(groups)[1]])) 
    max_phyto = round(log10(param$wMax_phyto), digits = 2)
    zoo_start =  which(round(log10(model$w), digits = 2) == start)
    zoo_finish = which(round(log10(model$w), digits = 2) == finish) 
    zoo_slope2 = lm(log10(tot_anim[zoo_start:zoo_finish])~log10(model$w[zoo_start:zoo_finish]))$
      coefficients[2]
    #zoo_slope2 =  round((log10(tot_zoo[zoo_finish]) - log10(tot_zoo[zoo_start]))/
    #                      (log10(model$w[zoo_finish]) - log10(model$w[zoo_start])), digits = 2)
    fish_slope2 = round(lm(log10(tot_fish[fish_start:fish_finish])~log10(model$w[fish_start:fish_finish]))$
                          coefficients[2], digits = 2)
    #fish_slope2 = round((log10(N_ave[dim(N_ave)[1],fish_finish]) - log10(N_ave[dim(N_ave)[1],fish_start]))/
    #                      (log10(model$w[fish_finish]) - log10(model$w[fish_start])), digits = 2)
    phyto_slope = enviro$b
    
    slopes[i,1] <- phyto_slope
    slopes[i,2] <- zoo_slope2
    slopes[i,3] <- fish_slope2
  }
  
slopes
}


### TROPHIC LEVEL
t_level = function(diets, enviro_data, group, en){
stor = matrix(0, nrow = dim(enviro_data)[1], ncol = 12)

### STARTING TROPHIC LEVELS
starter_tl = c(1, 2, 2, 2, 2,2,2, 2, 2, 2, 2, 2, 2)

for(i in 1:dim(enviro_data)[1]){
  start_tl = starter_tl
  diet = diets[[i]]
  phyto_diet = rowSums(diet[,c(1:3)])
  diet = cbind(phyto_diet, diet[,-c(1:3)])
  tester = round(100*(diet/rowSums(diet)))
  for(k in 1:100){
    for(j in 1:dim(tester)[1]){
      curr_tl =  sum(start_tl*tester[j,]/100)
      start_tl[j+1] = curr_tl + 1
    }
   }
  start_tl = start_tl[2:13]
  stor[i,] = start_tl
  
  }

namess = c("Flagellates", "Ciliates", "Larvaceans", "Omnivorous Copepods", "Carnivorous Copepods", "Euphausiids",
           "Chaetognaths", "Salps", "Jellyfish", "Fish (< 100g)", "Fish (100g-10kg)", "L Fish (>10kg)")

par(mfrow = c(3,3))
stor = stor[,c(1,2,3,4,6,8,5,7,9,10,11,12)]
namess = namess[c(1,2,3,4,6,8,5,7,9,10,11,12)]
for(i in 1:9){
if(en == "chlo"){
  if(i < 7){
plot(log10(enviro_data$chlo), stor[,i], ylab = "Trophic Level", main = namess[i],
     xlab = expression(paste("log"[10], "(Chlorophyll mg m"^-3, ")")), ylim = c(2,3))
  }
if(i > 6){
  plot(log10(enviro_data$chlo), stor[,i], ylab = "Trophic Level", main = namess[i],
       xlab = expression(paste("log"[10], "(Chlorophyll mg m"^-3, ")")), ylim = c(3.2,4))
}
    }
}
  
  par(mfrow = c(1,1))
 plot(log10(enviro_data$chlo), stor[,10], ylab = "Trophic Level", 
           xlab = expression(paste("log"[10], "(Chlorophyll mg m"^-3, ")")), ylim = c(3,6))
 points(log10(enviro_data$chlo), stor[,11], col = "red")
 points(log10(enviro_data$chlo), stor[,12], col = "blue")
 legend("topright", legend = c("Small", "Medium", "Large"), fill = c("black", "red", "blue"))
 
}


####### OLIGO VERSUS EUTRO
food_web <- function(N, enviro_data, oligo_chlo, eutro_chlo){

  ## BIOMASS OF THE GROUPS
  biom_mat <- matrix(0, nrow = dim(enviro_data)[1], ncol = 12)
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    zoo_abunds = rowSums(B_ave)
    zoo_props = zoo_abunds
    biom_mat[i,] <- zoo_props 
  }
  

  ## TROPHIC LEVEL OF GROUPS
  stor = matrix(0, nrow = dim(enviro_data)[1], ncol = 13)

  starter_tl = c(1, 2, 2, 2, 2, 3, 2, 3, 2, 3, 4, 4, 4)
  
  for(i in 1:dim(enviro_data)[1]){
    start_tl = starter_tl
    if(dim(diets[[i]])[1] == 13){
      tester = round(100*(diets[[i]]/rowSums(diets[[i]])))
      
      for(j in 1:dim(tester)[1]){
        curr_tl =  sum(start_tl*tester[j,]/100)
        start_tl[j+1] = curr_tl + 1
      }
      
    }
    start_tl = start_tl[2:14]
    stor[i,] = start_tl
  }
  
  namess = c("Phyto","Flagellates", "Ciliates", "Larvaceans", "Omnivorous Copepods", "Carnivorous Copepods", "Euphausiids",
             "Chaetognaths", "Salps", "Jellyfish", "S Fish", "M Fish", "L Fish")
  
  oligo = which(log10(enviro_data$chlo) <= oligo_chlo)
  eutro = which(log10(enviro_data$chlo) >= eutro_chlo)
  
  oligo_phyto =  mean((enviro_data[oligo,"chlo"]^0.89)*(10^1.79))
  eutro_phyto =  mean((enviro_data[eutro,"chlo"]^0.89)*(10^1.79))
    
  bio_tl = matrix(0, nrow = 13, ncol = 4)
  rownames(bio_tl) = namess
  bio_tl[,1] = c(oligo_phyto,colMeans(biom_mat[oligo,])*1000)
  bio_tl[,2] = c(1,colMeans(stor[oligo,c(1:12)]))
  bio_tl[,3] = c(eutro_phyto,colMeans(biom_mat[eutro,])*1000)
  bio_tl[,4] = c(1,colMeans(stor[eutro,c(1:12)]))

  colnames(bio_tl) = c("Oligo Biom", "Oligo TL", "Eutro Biom", "Eutro TL")

diet_store = matrix(0, nrow = 12, ncol = 13)
diet_array = round(array(as.numeric(unlist(diets)), dim = c(13,13,dim(enviro_data)[1])),5)*(1/enviro_data$dt[1])*1000

oligo = which(log10(enviro_data$chlo) <= oligo_chlo)
eutro = which(log10(enviro_data$chlo) >= eutro_chlo)

oligo_diet = apply(diet_array[,,oligo], c(1,2), mean)
eutro_diet = apply(diet_array[,,eutro], c(1,2), mean)

return(list("btl" = bio_tl, "od" = oligo_diet, "ed" = eutro_diet))
}


pre_pressure <- function(diets, enviro_data, groups){
  diet_pressure = matrix(0, nrow = dim(enviro_data)[1], ncol = 13)
  growth = matrix(0, nrow = dim(enviro_data)[1], ncol = 12)
  biomass_mat = matrix(0, nrow = dim(enviro_data)[1], ncol = 13)
  colnames(diet_pressure) = c("Phytoplankton",as.character(groups$species))
  temp_effect = exp(25.55 - 0.63/(8.62e-05*(273+round(enviro_data$sst, digits = 1))))
  
  for(i in 1:dim(enviro_data)[1]){
    curr_diet <- diets[[i]]
    assim_mat <- 0.25*cbind((0.15/groups$carbon), matrix(c(groups$carbon), nrow = 12, ncol = 12, byrow = TRUE)/
                  matrix(c(groups$carbon), nrow = 12, ncol = 12) )
    growth[i,] <- rowSums(curr_diet*assim_mat)/temp_effect[i]
    diet_pressure[i,] <- colSums(curr_diet)/temp_effect[i]
    N_ave = N[[i]]
    B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    bioms = rowSums(B_ave)
    biomass_mat[i,2:13] <- bioms
    biomass_mat[i,1] <- enviro_data$chlo[i]*60/(0.1*1000)
  }
  
  par(mfrow = c(4,3))
  for(i in 1:13){
    plot(log10(enviro_data$chlo),diet_pressure[,i]/biomass_mat[,i], main = colnames(diet_pressure)[i])
    abline(1,0, col = "red")
    abline(v = -0.5, col = "blue")
    }
  
  par(mfrow = c(2,2))
  for(i in c(3,4,6,8)){
    plot(log10(enviro_data$chlo),diet_pressure[,i+1]/growth[,i], main = colnames(diet_pressure)[i+1], ylab = "")
    abline(0,0, col = "red")
    abline(v = -0.5, col = "blue")
    abline(h = 1, col = "red")
  }
  title("Predation/Growth", outer = TRUE)
  
  par(mfrow = c(1,1))
  plot(log10(enviro_data$chlo), rowSums(growth)/rowSums(diet_pressure[,2:13]))
  
  who_eats_who = matrix(0, nrow = dim(enviro_data)[1], ncol = 12)
  
  for(i in 1:dim(enviro_data)[1]){
    curr_diet <- diets[[i]]
    diet_p <- round(curr_diet/matrix(colSums(curr_diet), nrow = 12, ncol = 13, byrow = TRUE),2)
    who_eats_who[i,] <- diet_p[,group]
  }
  
  num_plots = which(colSums(who_eats_who) > 0.01)
  num_rows = ceiling(length(num_plots)/2)
  
  par(mfrow = c(num_rows, 2))
  for(i in 1:length(num_plots)){
    plot(log10(enviro_data$chlo), who_eats_who[,num_plots[i]],
         main = groups$species[num_plots[i]], ylab = "Pred Press Prop")
    if(group > 1){
      title(groups$species[group-1], outer = TRUE)
    }
    if(group == 1){
      title("Phytoplankton", outer = TRUE)
    }
  }
  
}



### TROPHIC LEVEL
t_level_mod = function(diets, enviro_data, group, en){
  stor = matrix(0, nrow = dim(enviro_data)[1], ncol = 10)
  
  ### STARTING TROPHIC LEVELS
  starter_tl = c(1, 2,2,2,2, 2, 2, 2, 2, 2, 2)
  
  for(i in 1:dim(enviro_data)[1]){
    start_tl = starter_tl
    diet = diets[[i]]
    diet = diet[-c(3,8),-c(4,9)]
    tester = round(100*(diet/rowSums(diet)))
    for(k in 1:100){
      for(j in 1:dim(tester)[1]){
        curr_tl =  sum(start_tl*tester[j,]/100)
        start_tl[j+1] = curr_tl + 1
      }
    }
    start_tl = start_tl[2:11]
    stor[i,] = start_tl
    
  }
  
  namess = c("Flagellates", "Ciliates", "Omnivorous Copepods", "Carnivorous Copepods", "Euphausiids",
             "Chaetognaths", "Jellyfish", "Fish (< 100g)", "Fish (100g-10kg)", "L Fish (>10kg)")
  colnames(stor) = namess
  par(mfrow = c(3,3))
  for(i in 1:7){
    if(en == "chlo"){
      plot(log10(enviro_data$chlo), stor[,i], ylab = "Trophic Level", main = namess[i],
           xlab = expression(paste("log"[10], "(Chlorophyll mg m"^-3, ")")))
    }
  }
  
  par(mfrow = c(3,1))
  for(i in 8:10){
    if(en == "chlo"){
      plot(log10(enviro_data$chlo), stor[,i], ylab = "Trophic Level", main = namess[i], 
           xlab = expression(paste("log"[10], "(Chlorophyll mg m"^-3, ")")))
    }
    
  }
}


par(mfrow = c(3,1))
plot(log10(enviro_data$chlo), storr[,10])
plot(log10(enviro_data$chlo), storrr[,8])
plot(log10(enviro_data$chlo), storrr[,8]-storr[,10])
abline(0,0, col = "red", lwd = 2)

par(mfrow = c(3,1))
plot(log10(enviro_data$chlo), storr[,11])
plot(log10(enviro_data$chlo), storrr[,9])
plot(log10(enviro_data$chlo), storrr[,9]-storr[,11])
abline(0,0, col = "red", lwd = 2)

par(mfrow = c(3,1))
plot(log10(enviro_data$chlo), storr[,12])
plot(log10(enviro_data$chlo), storrr[,10])
plot(log10(enviro_data$chlo), storrr[,10]-storr[,12])
abline(0,0, col = "red", lwd = 2)





diet_pressure = matrix(0, nrow = dim(enviro_data)[1], ncol = 13)
growth = matrix(0, nrow = dim(enviro_data)[1], ncol = 12)

biomass_mat = matrix(0, nrow = dim(enviro_data)[1], ncol = 13)

phyto_max = enviro_data$phyto_max
phyto_max[phyto_max >= -7.2] <- -7.2

for(i in 1:dim(enviro_data)[1]){
  N_ave = N[[i]]
  B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
  biomass_mat[i,c(2:13)] = rowSums(B_ave)
  
  phyto_range = seq(-14.5, phyto_max[i], 0.1)
  inter = enviro_data$a[i]
  slope = enviro_data$b[i]
  
  phyto = sum((10^inter)*((10^phyto_range))^(slope + 1))
  biomass_mat[i,1] = phyto
}

plot(log10(enviro_data$chlo), 
     rowSums(biomass_mat[,c(2:10)])/biomass_mat[,1], ylab = "Z:P")

colnames(diet_pressure) = c("Phytoplankton",as.character(groups$species))
temp_effect = exp(25.55 - 0.63/(8.62e-05*(273+round(enviro_data$sst, digits = 1))))

for(i in 1:dim(enviro_data)[1]){
  curr_diet <- diets[[i]]
  assim_mat <- 0.25*cbind((0.15/groups$carbon), matrix(c(groups$carbon), nrow = 12, ncol = 12, byrow = TRUE)/
                            matrix(c(groups$carbon), nrow = 12, ncol = 12) )
  growth[i,] <- rowSums(curr_diet*assim_mat)/temp_effect[i]
  diet_pressure[i,] <- colSums(curr_diet)/temp_effect[i]
  N_ave = N[[i]]
  B_ave = N_ave*matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
  bioms = rowSums(B_ave)
  biomass_mat[i,2:13] <- bioms
  biomass_mat[i,1] <- enviro_data$chlo[i]*60/(0.1*1000)
}

par(mfrow = c(4,3))
for(i in 1:13){
  plot(log10(enviro_data$chlo),diet_pressure[,i]/biomass_mat[,i], main = colnames(diet_pressure)[i])
  abline(1,0, col = "red")
  abline(v = -0.5, col = "blue")
}


######### CALCULATE AVERAGE GROWTH RATES FOR EACH GROUP
growth_ave <- function(growth, N, enviro_data){
  
  growth_mat = matrix(0, nrow = dim(enviro_data)[1], ncol = 12)
  biom = matrix(0, nrow = dim(enviro_data)[1], ncol = 12)
  growth_array = matrix(0, nrow = 12, ncol = 168)
  
  for(i in 1:dim(enviro_data)[1]){
    N_ave = N[[i]]
    gg = growth[[i]]
    growth_array = growth_array + gg/(dim(enviro_data)[1])
    temp_effect = 2^((enviro_data$sst[i] - 15.6)/10)
    gg_daily = gg/(365*1)
    w = matrix(model$w, nrow = dim(N_ave)[1], ncol = dim(N_ave)[2], byrow = TRUE)
    gg_adj = gg_daily/w
    biom[i,] = rowMeans(w*N_ave)
    for(j in 1:12){
    growth_mat[i,j] <- sum(gg_adj[j,])/sum(c(gg_adj[j,] > 0))
    }
  }
  
  
  par(mfrow = c(4,2))
  for(i in 3:9){
  plot(log10(enviro_data$chlo), log10(growth_mat[,i]), main = Groups$species[i])  
  }
  
  grow_zoo = growth_mat[,3:12]
  biom_zoo = biom[,10:12]
  
  gg_table = data.frame("min" = apply(grow_zoo, 2, min), 
                        "max" = apply(grow_zoo, 2, max), 
                        "mean" = apply(grow_zoo, 2, mean),
                        row.names = Groups$species[3:12])
  
}

## 12 - oligo
## 82 - max eutro

gg = growth[[1]]/365

gg = growth_array/365



w_range = log10(w[1,c(1:138)])
gg_l = log10(gg[1:12,c(1:138)])
gg_l[is.infinite(gg_l)] = NA

coll = rainbow(12)
par(mfrow = c(1,1))
plot((w_range), (gg_l[1,]-w_range), 
     ylim = c(-3, 0),
     ylab = "log10(Growth Rate , day-1)", xlab = "log10(w)", type = "l", lwd = 2,
     col = coll[1], main = "Oligo")

for(i in 2:12){
  lines(w_range, c(gg_l[i,]-w_range), lwd = 2, col = coll[i])
}
abline(0,0, col = "black", lwd = 2, lty = 2)
leger = c(as.character(Groups$species)[1:9], "Fish")

legend("topright", legend = leger, lty = 1,col = coll[c(1:9,12)], bty = "n",
       cex = 0.7, lwd = 2)

## CONVERT TO CARBON SCALE
carbo = matrix(Groups$carbon, nrow = 12, ncol= 168)
carbo_gg = carbo*gg
shifter = round(log10(1/Groups$carbon),1)/0.1
shifterr = shifter - 8
carbo_gg2 = carbo_gg*0

for(i in 1:12){
  w_min = which(round(log10(w[1,]),2) == Groups$W0[i]) 
  w_max = which(round(log10(w[1,]),2) == Groups$Wmax[i]) 
  
  w_min_s = which(round(log10(w[1,]),2) == Groups$W0[i]) - shifterr[i]
  w_max_s = which(round(log10(w[1,]),2) == Groups$Wmax[i]) - shifterr[i]
  
  carbo_gg2[i,w_min_s:w_max_s] = carbo_gg[i,w_min_s:w_max_s]
}

carb_w = log10(w[1,]) - 0.8
carb_w = carb_w[1:168]

carbo_gg2_l = log10(carbo_gg2)
carbo_gg2_l[is.infinite(carbo_gg2_l)] = NA

coll = rainbow(12)
par(mfrow = c(1,1))
plot((carb_w), (carbo_gg2_l[1,]), 
     ylim = c(min(carbo_gg2_l, na.rm = TRUE), max(carbo_gg2_l, na.rm = TRUE)),
     ylab = "log10(Growth Rate, g C day-1)", xlab = "log10(g C)", type = "l", lwd = 2,
     col = coll[1], xlim = c(-11.5, 2), main = "Mean Growth Rates")

for(i in 2:12){
  lines(carb_w, carbo_gg2_l[i,], lwd = 2, col = coll[i])
}
abline(0,1, col = "black", lwd = 2, lty = 2)

growth_KH = read.csv("growth_KH.csv")
points(log10(growth_KH$g_C_ind), log10(growth_KH$g_C_d))

leger = c(as.character(Groups$species)[1:9], "Fish")

legend("bottomright", legend = leger, lty = 1,col = coll[c(1:9,12)], bty = "n",
       cex = 0.7, lwd = 2)


