## Function to calculate slope intercept and maximum size of phytoplankton spectrum, for zooplankton
## resolved size spectrum model (Heneghan et al. in prep).

# Last updated March 2020

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
  df$phyto_max <- min(-7, df$phyto_max)

  return(df)
}

## OLD CODE USING HIRATA ET AL. (2011)

#
# ## Function to calculate slope and intercept of phytoplankton spectrum,
# # using empirical model from Hirata et al., (2011)
#
# phyto_info = function(chlo){ # chlo is chlorophyll concentration in mg m^-3
#
#   ## Calculate pico, nano, micro phytoplankton proportions
#   pico <- -(0.1529 + exp(1.0306*log10(chlo) - 1.5576))^(-1) - 1.8597*log10(chlo) + 2.9954
#   micro <- (0.9117 + exp(-2.7330*log10(chlo) + 0.4003))^(-1)
#   nano <- 1 - micro - pico
#   pico[pico > 1]  = micro[micro > 1] = 1
#   pico[pico < 0] = 0
#   nano[nano < 0] = 0
#   micro[micro < 0] = 0
#
#   ## Convert total chlorophyll to g m^-3 total wet weight - biomass
#   ## Allocate total chlorophyll to the three size classes
#   c_chl <- ((chlo^0.89)*(10^1.79))/chlo # chlo:carbon ratio (0.02 Chl:C)
#   tot_biom <- c_chl*chlo/1000 # (convert to grams carbon)
#
#   pico_biom <- pico*tot_biom
#   nano_biom <- nano*tot_biom
#   micro_biom <- micro*tot_biom
#
#   ## Convert from biomass density to abundance density
#   ## Convert from grams carbon to grams wet weight
#   w_max = 0.1*round((-8.4 + 4*micro)/0.1) # Maximum size depends on the proportion of micro
#   w_max = min(-5.5, w_max)
#   pico_range = seq(-14.5, -11.5, 0.1) # (0.2-2um)
#   nano_range = seq(-11.4, -8.5, 0.1) # (2-20um)
#   micro_range = seq(-8.4, w_max, 0.1) # (20-100um)
#
#   # Wet weight to carbon ratio
#   ww_carb_p <- 10 #1/((0.216*(1e12*(10^pico_range))^0.939/1e12)/(10^pico_range))
#   ww_carb_n <- 10 #1/((0.216*(1e12*(10^nano_range))^0.939/1e12)/(10^nano_range))
#   ww_carb_m <- 10 #1/((0.216*(1e12*(10^micro_range))^0.939/1e12)/(10^micro_range))
#
#   pico_biom <- (pico_biom/length(pico_range))*ww_carb_p
#   nano_biom <- (nano_biom/length(nano_range))*ww_carb_n
#   micro_biom <- (micro_biom/length(micro_range))*ww_carb_m
#
#   # Find abundance in each size range by assuming equal biomass in 0.1 size bins, then dividing by
#   # average individual size in the size bing
#   pico_abund = pico_biom/(10^pico_range)
#   nano_abund = nano_biom/(10^nano_range)
#   micro_abund = micro_biom/(10^micro_range)
#
#   na.action = na.omit
#
#   y = c(log10(pico_abund), log10(nano_abund), log10(micro_abund))
#   y[is.infinite(y)] = NA
#   x = c((pico_range), nano_range, micro_range)
#
#   log_int = summary(lm(y~x))$coefficients[1]
#   slope = summary(lm(y~x))$coefficients[2]
#   #plot(x,y)
#   #abline(log_int, slope)
#   pico_biomm <- sum((10^log_int)*((10^pico_range))^(slope + 1))
#   nano_biomm <- sum((10^log_int)*((10^nano_range))^(slope + 1))
#   micro_biomm <- sum((10^log_int)*((10^micro_range))^(slope + 1))
#   model_biom <- sum(micro_biomm + nano_biomm + pico_biomm)
#   actual_biom <- sum(length(micro_range)*(micro_biom) + length(nano_range)*sum(nano_biom) +
#                        length(pico_range)*sum(pico_biom))
#   return(list(log_int, slope, w_max, model_biom, actual_biom))
# }
