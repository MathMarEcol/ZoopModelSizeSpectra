### Author: Ryan Heneghan
### Date: May 2019
### This script creates figures comparing output from the size spectrum model, with
### output from the statistical models for different zooplankton groups.

rm(list = ls())
setwd("~/Dropbox/EcosystemEfficiency/NumericalModel/ModelFigures/Figure Code")

library(ggplot2)
library(egg)
library(maptools)
data(wrld_simpl)

############# ############# ############# ############# ############# 
############# ############# ############# ############# ############# 
############# ABUNDANCE EMPIRICAL VERSUS MODEL PLOTS

# Lat ranges of the statistical model output (where 85% of the data is)
lat_ranges = data.frame("min" = c(-40, -48, -43, -50, -42, -38, -40), "max" = c(40, 48, 43, 50, 42, 38, 40))

th_abunds = read.csv("th_abunds_intt.csv") ## Import abundances from size spectrum model (5x5 degree grid square resolution)
st_abunds = read.csv("5deg_year_ave_data.csv") ## Import abundances from statistical models (aggregated to annual, 5x5 degree grid squares)
enviro_data = read.csv("enviro_5d_new.csv") ## Import environmental data
names(enviro_data)[names(enviro_data) == "x"] <- "lon" # Rename x variable to lon for longitude
names(enviro_data)[names(enviro_data) == "y"] <- "lat" # Rename y variable to lat for latitude

# remove land grid squares
remove_these <- which(is.na(enviro_data$chlo) == TRUE | is.na(enviro_data$sst) == TRUE) #identify land grid squares - where is sst and chlo NA?
enviro_data = enviro_data[-remove_these,]

# Remove these from st_abunds (there may be some stray coastal cells that exist after the statistical model aggregation, this will get rid of them)
st_abunds = st_abunds[-remove_these,]

st_abunds = st_abunds[,-c(1,2)] # Remove lat and lon columns from statistical model output
st_abunds <- 10^(st_abunds) # Raise statistical output to power of 10 (if logged)

## Reorder chaetognaths and krill
th_abunds <- th_abunds[,c(1,2,3,5,4,6,7)]
st_abunds <- st_abunds[,c(1,2,3,5,4,6,7)]

st_props <- st_abunds/rowSums(st_abunds) # Calculate proportions from statistical model output
th_props <- th_abunds/rowSums(th_abunds) # Calculate proportions from size spectrum model output

############# ############# ############# ############# ############# 
#################  PLOT 1
library(maptools)
data(wrld_simpl)

plot_list = list()
l_brack = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)", "j)", "k)", "l)")
plot_names = c("Larvaceans", "Omnivorous\nCopepods", "Carnivorous\nCopepods", "Chaetognaths")

for(x in 1:4){ # This for loop does the first four zooplankton groups (larvs, occops, ccops, chaets)
  ## Set lat limits, define lat and long
  lat_min = lat_ranges[x, "min"]
  lat_max = lat_ranges[x, "max"]
  lat_lons = c(enviro_data$lat < lat_max & enviro_data$lat > lat_min)
  lat_long = enviro_data[lat_lons, c("lat", "lon")]
  
  ## Absolute abundances
  st_curr = log10(st_abunds[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x])
  th_curr = log10(th_abunds[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x])
  
  ## Abundance proportions  
  #st_curr = st_props[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  #th_curr = th_props[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  
  ## Create data frame for ggplot
  st_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = st_curr)
  th_dat = data.frame("lat" = lat_long$lat, "lon" = lat_long$lon, "zoo" = th_curr)
  
  ## Empirical map
  st_plot = ggplot(data = st_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster()  + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_text(size =18, margin = unit(c(0,0.2,0,0.2), unit = "cm")),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(-.8,0,0.0,0), "cm"),
                       legend.text = element_text(size = 14),
                       legend.box.margin = margin(0,0,0,-10),
                       panel.grid = element_blank(),
                       plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab(plot_names[x]) + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[x + (x-1)*2])
  
  ## Size spectrum map
  th_plot = ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(-.8,0,0.0,0), "cm"),
                       plot.subtitle=element_text(size=16, color="black"),
                       legend.text = element_text(size = 14),
                       panel.grid = element_blank(),
                       legend.box.margin = margin(0,0,0,-10)) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab("") + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[(x+1) + (x-1)*2])
  
  ## Calculate correlations
  curr_data = data.frame("x" = st_curr, "y" = th_curr)
  curr_cor = round(cor(st_curr, th_curr, use = "pairwise.complete.obs"), digits = 2)
  
  ## Set legend position for correlation plot
  curr_legend = bquote(rho ~ " = " ~.(curr_cor))
  curr_x_pos = min(st_curr, na.rm = TRUE) + 0.75*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  curr_y_pos = min(th_curr, na.rm = TRUE) + 0.15*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  #label_data = data.frame("x" = curr_x_pos, "y" = curr_y_pos, "curr_legend" = expression(alpha))
  
  coorder = 0.6*(max(st_curr, na.rm = TRUE) - min(st_curr, na.rm = TRUE))/(max(th_curr, na.rm = TRUE) - min(th_curr, na.rm = TRUE))
  axmin = curr_x_pos - 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  axmax = curr_x_pos + 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  aymin = curr_y_pos - 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  aymax = curr_y_pos + 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  tt=0.01
  l = paste("r ==", curr_cor)
  
  ## Correlation plot
  corr_plot = ggplot(curr_data, aes(x, y)) + geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    xlab("Empirical") + ylab("Theoretical") + 
    theme_classic() + theme(plot.margin = unit(c(-.5,0.5,0.0,0), "cm"), plot.subtitle=element_text(size=16, color="black"),
                            axis.text = element_text(size = 14), axis.title = element_text(size = 14)) + coord_fixed(coorder) +
    scale_x_continuous(limits = c(min(st_curr, na.rm = TRUE), max(st_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    scale_y_continuous(limits = c(min(th_curr, na.rm = TRUE), max(th_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    ggtitle(label = "",subtitle = l_brack[(x+2) + (x-1)*2]) + 
    annotate("rect", xmin = axmin, xmax = axmax,  ymin = aymin, ymax = aymax, alpha = "1", fill = "gray", color = "black") +
    annotate("text", x = c(curr_x_pos), y = c(curr_y_pos), parse = TRUE, label = l, size = 5)
  
  if(x == 1){ ## Set column titles
    st_plot = st_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.4,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop('Empirical', paste(log[10],'(# ', m^-3,')'))), subtitle = "a)")
    
    th_plot = th_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.4,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop('SS Model', paste(log[10],'(# ', m^-3,')'))), subtitle = "b)")
    
    corr_plot = corr_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(0.2,0,0,0), unit = "cm")),
                                  plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop("Empirical v SS Model", "Correlation (r)")), subtitle = "c)")
    
  }
  
  ## Save plots to lists
  plot_list[[x + (x-1)*2]] <- st_plot
  plot_list[[(x+1) + (x-1)*2]] <- th_plot
  plot_list[[(x+2) + (x-1)*2]] <- corr_plot
  
}

## Plot figures
ggsave(filename = "Figure2_standabs.png", plot = ggarrange(plots = plot_list, cliop = "off",labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)", "j)", "k)", "l)"), nrow = 4), width = 14, height = 12.4)

############# ############# ############# 
##############  PART 2
plot_list = list()
l_brack = c("m)", "n)", "o)", "p)", "q)", "r)", "s)", "t)", "u)")
plot_names = c("Euphausiids", "Salps", "Jellyfish")

for(x in 5:7){ # This for loop does the last three zooplankton groups (euphs, salps, jellys)
  ## Set lat limits, define lat and long
  lat_min = lat_ranges[x, "min"]
  lat_max = lat_ranges[x, "max"]
  lat_lons = c(enviro_data$lat < lat_max & enviro_data$lat > lat_min)
  lat_long = enviro_data[lat_lons, c("lat", "lon")]
  st_curr = log10(st_abunds[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x])
  th_curr = log10(th_abunds[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x])
  
  ## Absolute abundances
  st_curr = log10(st_abunds[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x])
  th_curr = log10(th_abunds[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x])
  
  ## Abundance proportions  
  #st_curr = st_props[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  #th_curr = th_props[c(enviro_data$lat < lat_max & enviro_data$lat > lat_min), x]
  
  ## Create data frame for ggplot
  curr_xer = x - 4 # curr_xer for calling plot names
  
  # Empirical plot
  st_plot = ggplot(data = st_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster()  + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_text(size =18, margin = unit(c(0,0.2,0,0), unit = "cm")),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(-.8,0,0.0,0), "cm"),
                       legend.text = element_text(size = 14),
                       legend.box.margin = margin(0,0,0,-10),
                       panel.grid = element_blank(),
                       plot.subtitle=element_text(size=16, color="black")) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab(plot_names[curr_xer]) + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[curr_xer + (curr_xer-1)*2])
  
  # Size spectrum plot
  th_plot = ggplot(data = th_dat, aes(x = lon, y = lat, fill = zoo))  + 
    geom_raster() + scale_fill_gradientn(colours=c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), guide="colorbar", na.value = "white") + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey') + 
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       plot.margin = unit(c(-.8,0,0.0,0), "cm"),
                       plot.subtitle=element_text(size=16, color="black"),
                       legend.text = element_text(size = 14),
                       panel.grid = element_blank(),
                       legend.box.margin = margin(0,0,0,-10)) + coord_fixed(ratio = 1.) + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group=group), fill='grey', colour = "black") + 
    xlab("") + ylab("") + labs(fill = "") +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = "",subtitle = l_brack[(curr_xer+1) + (curr_xer-1)*2])
  
  ## Correlation plot
  curr_data = data.frame("x" = st_curr, "y" = th_curr)
  curr_cor = round(cor(st_curr, th_curr, use = "pairwise.complete.obs"), digits = 2)
  
  ## Set position and size of legend for correlation plot
  curr_legend = bquote(rho ~ " = " ~.(curr_cor))
  curr_x_pos = min(st_curr, na.rm = TRUE) + 0.80*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  curr_y_pos = min(th_curr, na.rm = TRUE) + 0.15*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  #label_data = data.frame("x" = curr_x_pos, "y" = curr_y_pos, "curr_legend" = expression(alpha))
  
  coorder = 0.6*(max(st_curr, na.rm = TRUE) - min(st_curr, na.rm = TRUE))/(max(th_curr, na.rm = TRUE) - min(th_curr, na.rm = TRUE))
  axmin = curr_x_pos - 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  axmax = curr_x_pos + 0.15*(max(st_curr, na.rm = TRUE)-min(st_curr,na.rm = TRUE))
  aymin = curr_y_pos - 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  aymax = curr_y_pos + 0.09*(max(th_curr, na.rm = TRUE)-min(th_curr,na.rm = TRUE))
  tt=0.01
  l = paste("r ==", curr_cor)
  
  ## Plot correlations
  corr_plot = ggplot(curr_data, aes(x, y)) + geom_point() + 
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    xlab("Empirical") + ylab("Theoretical") + 
    theme_classic() + theme(plot.margin = unit(c(-.8,0,0.0,0), "cm"), plot.subtitle=element_text(size=16, color="black"),
                            axis.text = element_text(size = 14), axis.title = element_text(size = 14)) + coord_fixed(coorder) +
    scale_x_continuous(limits = c(min(st_curr, na.rm = TRUE), max(st_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    scale_y_continuous(limits = c(min(th_curr, na.rm = TRUE), max(th_curr, na.rm = TRUE)), expand = c(0.01,0.01)) +
    ggtitle(label = "",subtitle = l_brack[(curr_xer+2) + (curr_xer-1)*2]) + 
    annotate("rect", xmin = axmin, xmax = axmax,  ymin = aymin, ymax = aymax, alpha = "1", fill = "gray", color = "black") +
    annotate("text", x = c(curr_x_pos), y = c(curr_y_pos), parse = TRUE, label = l, size = 5)
  
  if(curr_xer == 1){ # Set column names
    st_plot = st_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(1,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop('Empirical Model', paste(log[10],'(# ', m^-3,')'))), subtitle = "m)")
    
    th_plot = th_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(1,0,0,0), unit = "cm")),
                              plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop('SS Model', paste(log[10],'(# ', m^-3,')'))), subtitle = "n)")
    
    corr_plot = corr_plot + theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = unit(c(1,0,0,0), unit = "cm")),
                                  plot.subtitle=element_text(size=14, color="black")) + 
      ggtitle(expression(atop("Empirical v SS Model", "Correlation (r)")), subtitle = "o)")
    
  }
  
  # Save figures
  plot_list[[curr_xer + (curr_xer-1)*2]] <- st_plot
  plot_list[[(curr_xer+1) + (curr_xer-1)*2]] <- th_plot
  plot_list[[(curr_xer+2) + (curr_xer-1)*2]] <- corr_plot
  
}

# Plot figures
ggsave(filename = "figure2_standabs2.png", plot = ggarrange(plots = plot_list, nrow = 3), width = 14, height = 9)

