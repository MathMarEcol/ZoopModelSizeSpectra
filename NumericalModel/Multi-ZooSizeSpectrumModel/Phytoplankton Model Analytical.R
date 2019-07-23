## Function to calculate slope and intercept of phytoplankton spectrum, 
# using empirical model from Hirata et al., (2011)

phyto_info_2 = function(chlo){ # chlo is chlorophyll concentration in mg m^-3
  
  ## Calculate pico, nano, micro phytoplankton proportions
  ## BREWIN ET AL., 2015
  pico <- (0.13*(1-exp(-0.8/0.13*chlo)))/chlo
  nano <- (0.77*(1-exp(-0.94/0.77*chlo)))/chlo - pico
  micro <- (chlo - 0.77*(1-exp(-0.94/0.77*chlo)))/chlo
  
  ## Calculate maximum size
  w_max = 0.1*round((-8.4 + 3*micro)/0.1) # Maximum size depends on the proportion of micro
  w_max = min(-5.4, w_max)
  
  ## Convert total chlorophyll to g m^-3 total wet weight - biomass
  ## Allocate total chlorophyll to the three size classes
  c_chl <- ((chlo^0.89)*(10^1.79))/chlo # chlo:carbon ratio (0.02 Chl:C)
  tot_biom_c <- c_chl*chlo/1000 # (convert to grams carbon)
  tot_biom <- tot_biom_c*(1/0.1) # convert to grams wet weight, assuming 0.1 C:ww
  
  # Break up total biom into pico, nano and micro
  pico_biom <- pico*tot_biom
  nano_biom <- nano*tot_biom
  micro_biom <- micro*tot_biom

  ## Find abundances at boundaries of pico, nano size ranges, by analytically
  ## solving integral of N = aw^b
  
  pp = 2*pico_biom/(10^-11.5 - 10^-14.5)
  nn = 2*nano_biom/(10^-8.5 - 10^-11.5)

  N_1 = nn*pp/(nn + pp) # Abundance at 2um
  N_0 = pp - N_1 # Abundance at 0.2um
  N_2 = nn - N_1 # Abundance at 20um
  N_3 = 10^((w_max - -8.5)/(-8.5 -- 11.5)*(log10(N_2)-log10(N_1)) + log10(N_2))
  
  w_0 = -14.5
  w_1 = -11.5
  w_2 = -8.5
    
  b = (log10(pico_biom) - log10(nano_biom) - w_1 + w_2)/(w_1 - w_2)  
  a = pico_biom*(b+1)/((10^(w_1))^(b+1) - (10^(w_0))^(b+1))
  
  ## CHECK BIOMASS
  actual_biomass = pico_biom + nano_biom + micro_biom
  under_curve = a/(b+1)*((10^w_max)^(b+1) - (10^-14.5)^(b+1))
  pico_curve = a/(b+1)*((10^-11.5)^(b+1) - (10^-14.5)^(b+1))
  nano_curve = a/(b+1)*((10^-8.5)^(b+1) - (10^-11.5)^(b+1))
  micro_curve = a/(b+1)*((10^w_max)^(b+1) - (10^-8.5)^(b+1))
  valid_table = data.frame("under_curve" = c(pico_curve, nano_curve, micro_curve, under_curve),
                           "actual" = c(pico_biom, nano_biom, micro_biom, actual_biomass),
                           "curve/actual" = c(pico_curve, nano_curve, micro_curve, under_curve)/
                                         c(pico_biom, nano_biom, micro_biom, actual_biomass),
                           row.names = c("pico", "nano", "micro", "total"))
  return(list(a, b,under_curve, w_max))
}

