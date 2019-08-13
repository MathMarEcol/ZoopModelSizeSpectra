model_zoo <- function(curr_gam, BATH, sst_c, chlo_c, GEAR_MESH, adder){
  
  zoo_gam = curr_gam
  BATHY = BATH
  sst_clim = sst_c
  chlo_clim = chlo_c
  
  ### If we've asked for the maximum abundance gear_mesh level, find it
  if(GEAR_MESH == 'max'){
    GEAR_MESHER = names(zoo_gam$coefficients[grep('GEAR_MESH', names(zoo_gam$coefficients))])[which.max(zoo_gam$coefficients[grep('GEAR_MESH',names(zoo_gam$coefficients))])]
    GEAR_MESHER = unlist(strsplit(GEAR_MESHER, 'GEAR_MESH'))[2]
    adder = 0
  }
  
  ### Calculate mean gear_mesh level, then we feed the function the first gear_mesh level, and add the difference 'adder'
  ### to the output to get the abundances averaging across all gear_mesh levels
  if(GEAR_MESH == 'mean'){
    mean_gear_mesh = mean(zoo_gam$coefficients[grep('GEAR_MESH', names(zoo_gam$coefficients))], na.rm = TRUE)
    GEAR_MESHER = names(zoo_gam$coefficients)[grep('GEAR_MESH', names(zoo_gam$coefficients))[1]]
    GEAR_MESHER = unlist(strsplit(GEAR_MESHER, 'GEAR_MESH'))[2]
    adder = mean_gear_mesh - zoo_gam$coefficients[grep('GEAR_MESH', names(zoo_gam$coefficients))][1]
  }
  
  GEAR_MESH = GEAR_MESHER
  
  new_data <- data.frame(lon = BATHY$lon, lat = BATHY$lat, BATHY = abs(BATHY$bathy), 
                         GEAR_MESH = GEAR_MESH)
  
  month <- c("January", "February", "March", "April", "May", "June", "July", "August",
             "September", "October", "November", "December")
  
  results <- data.frame(lon = BATHY$lon, lat = BATHY$lat)
  results[,month] <- NA
  
  south_hemi <- c(BATHY$lat < 0)
  north_hemi <- c(BATHY$lat >= 0)
  
  if(sum(sst_clim[,c(3:14)], na.rm = TRUE) > 0){ # If sst and day of year are being used
  for(i in 1:12){
    n_day <- 15 + (i-1)*30 # Day of the year (approx. middle of month i) NORTH HEMISPHERE
    s_day <- (n_day + 183) %% 365 # SOUTH HEMPISPHERE DAY OF YEAR

    curr_month <- month[i]
    
    year_day <- c(rep(n_day, sum(north_hemi)), rep(s_day, sum(south_hemi)))
    month_data <- data.frame(SST = sst_clim[,curr_month], log10CHLO = log10(chlo_clim[,curr_month]), day_of_year = year_day)
    mon_new_data <- cbind(new_data, month_data)
    colnames(mon_new_data) <- c("lon", "lat", "BATHY", "GEAR_MESH", "SST","log10CHLO", "day_of_year")
    results[,curr_month] <- as.numeric(predict.gam(zoo_gam, mon_new_data)) + adder
    
  }
    results$year_ave <- rowMeans(results[,month], na.rm = TRUE)
    
  }
  
  if(sum(sst_clim[,c(3:14)], na.rm = TRUE) == 0){ # If we are holding sst and day of year constant
    results <- data.frame(lon = BATHY$lon, lat = BATHY$lat)
    results[,'year_ave'] <- NA
    
    year_ave_chlo = rowMeans(chlo_clim[,c(3:14)], na.rm = TRUE)
    
    year_day <- c(rep(n_day, sum(north_hemi)), rep(s_day, sum(south_hemi)))
    month_data <- data.frame(SST = sst_clim[,3], log10CHLO = log10(year_ave_chlo), day_of_year = 1)
    mon_new_data <- cbind(new_data, month_data)
    colnames(mon_new_data) <- c("lon", "lat", "BATHY", "GEAR_MESH", "SST","log10CHLO", "day_of_year")
    results[,'year_ave'] <- as.numeric(predict.gam(zoo_gam, mon_new_data)) + adder
    
  }
  return(results)
}


plot_all_abunds <- function(abund_data, lat_ranges, file_name1, file_name2){
  library("ggplot2")
  library("egg")
  library("maptools")
  data(wrld_simpl)
  
  plot_list = list()
  l_brack = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)", "j)", "k)", "l)", "m)", "n)")
  plot_names = colnames(abund_data)
  prop_data = abund_data[,-c(1,2)]/rowSums(abund_data[,-c(1,2)])
  
  for(x in 1:c(dim(abund_data)[2] - 2)){
    lat_min = lat_ranges[x, "min"]
    lat_max = lat_ranges[x, "max"]
    lat_lons = c(abund_data$lat < lat_max & abund_data$lat > lat_min)
    lat_long = abund_data[lat_lons, c("lat", "lon")]
    ab_curr = abund_data[lat_lons, c(x+2)]
    prop_curr = prop_data[lat_lons, x]
    
    ab_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = log10(ab_curr))
    pr_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = prop_curr)
    ab_plot = ggplot(data = ab_dat, aes(x = lon, y = lat, fill = zoo))  + 
      geom_raster()  + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
      theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                         axis.title.y = element_text(size =18, margin = unit(c(0,0.5,0,0), unit = "cm")),
                         axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                         panel.border = element_rect(colour = "black"),
                         plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
                         legend.text = element_text(size = 14),
                         plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
      geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
      xlab("") + ylab(plot_names[x+2]) + labs(fill = "") +
      scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
      scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
      ggtitle(label = "",subtitle = l_brack[(x-1)*2 + 1])
    
    max_colkey = max(prop_curr)
    pr_plot = ggplot(data = pr_dat, aes(x = lon, y = lat, fill = zoo))  + 
      geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white",limits = c(0, max_colkey)) + 
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
      ggtitle(label = "",subtitle = l_brack[(x-1)*2 + 2])
    
    if(x == 1){
      ab_plot = ab_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                                plot.subtitle=element_text(size=14, color="black")) + 
        ggtitle(expression(atop('Empirical Abundance', paste(log[10],'(# ', m^-3,')'))), subtitle = "a)")
      
      pr_plot = pr_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                                plot.subtitle=element_text(size=14, color="black")) + 
        ggtitle("Proportion of \n Abundance", subtitle = "b)")
    }
    
    if(x == 5){
      ab_plot = ab_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                                plot.subtitle=element_text(size=14, color="black")) + 
        ggtitle(expression(atop('Empirical Abundance', paste(log[10],'(# ', m^-3,')'))), subtitle = "i)")
      
      pr_plot = pr_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                                plot.subtitle=element_text(size=14, color="black")) + 
        ggtitle("Proportion of \n Abundance", subtitle = "j)")
    }
    
    plot_list[[(x-1)*2 + 1]] <- ab_plot
    plot_list[[(x-1)*2 + 2]] <- pr_plot
    
  }
  return(plot_list)

}
  
