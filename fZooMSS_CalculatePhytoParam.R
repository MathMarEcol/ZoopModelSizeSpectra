## Function to calculate slope intercept and maximum size of phytoplankton spectrum, for zooplankton
## resolved size spectrum model (Heneghan et al. in prep).

# Last updated 17th September 2020

fZooMSS_CalculatePhytoParam = function(df){ # chlo is chlorophyll concentration in mg m^-3

  ## Calculate pico, nano, micro phytoplankton proportions of total chlorophyll
  ## BREWIN ET AL., 2015
  pico <- (0.13*(1-exp(-0.8/0.13*df$chlo)))/df$chlo
  nano <- (0.77*(1-exp(-0.94/0.77*df$chlo)))/df$chlo - pico
  micro <- (df$chlo - 0.77*(1-exp(-0.94/0.77*df$chlo)))/df$chlo

  ## Convert total chlorophyll to g m^-3 total wet weight - biomass
  ## Allocate total chlorophyll to the three size classes
  c_chl <- ((df$chlo^0.89)*(10^1.79))/df$chlo # chlo:carbon ratio, from Mara??on et al. 2014
  tot_biom_c <- c_chl*df$chlo/1000 # (convert to grams carbon)
  tot_biom <- tot_biom_c*(1/0.1) # convert to grams wet weight, assuming 0.1 C:ww

  # Break up total biom into pico, nano and micro
  df$pico_biom <- pico*tot_biom
  df$nano_biom <- nano*tot_biom
  df$micro_biom <- micro*tot_biom

  ## Find abundances at boundaries of pico, nano size ranges, by analytically
  ## solving integral of N = aw^b

  w_0 <- -14.5 # log minimum size of picophytoplankton
  w_1 <- -11.5 # log minimum size of nanophytoplankton (max size of pico also)
  w_2 <- -8.5 # log minimum size of macrophytoplankton (max size of nano also)

  df$phyto_slope <- (log10(df$pico_biom) - log10(df$nano_biom) - w_1 + w_2)/(w_1 - w_2)  # Calculate slope
  df$phyto_int <- log10(df$pico_biom*(df$phyto_slope+1)/((10^(w_1))^(df$phyto_slope+1) - (10^(w_0))^(df$phyto_slope+1))) # Calculate intercept

  ## Calculate maximum size
  df$phyto_max <- 0.1*round((-8.4 + 2*micro)/0.1) # Maximum size depends on the proportion of micro
  max_phyto <- rep(-7, length(df$chlo)) # Set -7 to be the max possible size for phyto
  df$phyto_max <- pmin(max_phyto, df$phyto_max)

  return(df)
}
