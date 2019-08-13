plot_zoo <- function(zoo_data, name, save_plots){
  
data = zoo_data
zoo_name = as.character(name)
min_val = min(data[data$TOT_ABUND > 0, 'TOT_ABUND'], na.rm =TRUE)/2
no_obs = dim(data)[1]

if(save_plots == TRUE){
  # Plot location of data for data
  print(paste("Working on global map of", zoo_name, "observations"))
  graphics.off()
  png(filename= paste("ObsMap_", zoo_name, ".png", sep = ""), units="in",  width=7.5, height=6, res=200)
  library(ggplot2)
  library(maptools)
  data(wrld_simpl)
  st_plot <- ggplot(data = data, aes(x = LONGITDE, y = LATITUDE, colour = 'black'))  + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group = group), fill='grey', colour = "grey")+
    geom_point(colour = "black", size = 0.1)+
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       panel.grid = element_blank(),
                       plot.margin = unit(c(0,0.2,0.,0.2), "cm"),
                       legend.position = "none",
                       plot.subtitle=element_text(size=12, color="black")) + coord_fixed(ratio = 1.)+ 
    xlab("") + ylab("")  +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = paste(zoo_name, " observations", " (n = ", no_obs, ")", sep = ""))
  print(st_plot)
  dev.off()
    
print("Working on first set of histograms")
## Histograms of year, month, time of day, lat, chlo, bathymetry
png(filename = paste("Hist1_", zoo_name, ".png", sep = ""), units = "in", width = 7.5, height = 7, res = 200)
par(mfrow = c(2,2))
hist(data$YEAR,breaks = 50, ylab = "# Observations", xlab = "Year", main = "Years")
hist(as.numeric(data$MON), ylab = "# Observations", xlab = "months", main = "months (standardised to North Hemis)")
#if(sum(!is.na(as.numeric(data$TIMEloc)) > 0)){
#hist(as.numeric(data$TIMEloc), ylab = "# Observations", xlab = "Local Time of Day", main = "Time of Observations")
#}
hist(data$LATITUDE, ylab = "# Observations", xlab = "Latitude", main = "Latitude")
dev.off()

print("Working on second set of histograms")
png(filename =  paste("Hist2_", zoo_name, ".png", sep = ""), units = "in", width = 7.5, height = 7, res = 200)
par(mfrow = c(2,2))
hist(log10(data$CHLO), ylab = "# Observations", xlab = expression(paste("log"[10], "(Chlo-a, mg m"^-3, ")")), 
     main = "Chlorophyll")
hist(data$SST, ylab = "# Observations", xlab = "SST", main = "SST")    
hist(data$BATHY, ylab = "# Observations", xlab = "Bathymetry", main = "Bathymetry")
hist(log10(data$TOT_ABUND+min_val), ylab = "# Observations", xlab = expression(paste("log"[10], "(# datam"^-3, ")")), 
     main = "# m^-3")
dev.off()

print("Working on first set of abundance plots")
png(filename =  paste("Abund_Plots1_", zoo_name, ".png", sep = ""), units = "in", width = 6, height = 6, res = 200)
par(mfrow = c(2,2))
plot(data$BATHY, log10(data$TOT_ABUND+min_val), xlab = "Depth (metres)", 
     ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
     main = "Bathymetry", cex.lab = 1.3)
plot(data$SST, log10(data$TOT_ABUND+min_val), xlab = "SST", 
     ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
     main = "SST", cex.lab = 1.3)
plot(log10(data$CHLO), log10(data$TOT_ABUND+min_val), 
     xlab =expression(paste("log"[10], "(Chlo-a, mg m"^-3, ")")) , 
     ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
     main = "Chlorophyll", cex.lab = 1.3)
plot(data$MESH, log10(data$TOT_ABUND+min_val), xlab = "Mesh Size (um)", 
     ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
     main = "Mesh Size", cex.lab = 1.3)
dev.off()

print("Working on second set of abundance plots")
png(filename =  paste("Abund_Plots2_", zoo_name, ".png", sep = ""), units = "in", width = 6, height = 6, res = 200)
par(mfrow = c(3,2))
plot(data$MON, log10(data$TOT_ABUND+min_val), xlab = "month", 
     ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
     main = "month (North. Hemis.)", cex.lab = 1.3)
if(sum(!is.na(as.numeric(data$TIMEloc))) > 0){
  plot(as.numeric(data$TIMEloc), log10(data$TOT_ABUND+min_val), 
       xlab = "Local Hour" , 
       ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
       main = "Time of Day", cex.lab = 1.3)
}
plot(data$MID_Z, log10(data$TOT_ABUND+min_val),
     xlab = "Max. Tow Depth",
     ylab = expression(paste("log"[10], "(# m"^-3, ")")),
     main = "Sample Depth")
plot(data$SST, log10(data$CHLO),xlab = "SST", 
     ylab = expression(paste("log"[10], "(Chlo-a, mg m"^-3, ")")) , 
     main = "CHLO & SST", cex.lab = 1.3)
plot(data$MON, data$SST, 
     xlab = "month (North. Hemis.)" , 
     ylab = "SST", 
     main = "SST and month", cex.lab = 1.3)
dev.off()
}

if(save_plots == FALSE){
    
  print("Working on first set of histograms")
  ## Histograms of year, month, time of day, lat, chlo, bathymetry
  par(mfrow = c(2,2))
  hist(data$YEAR,breaks = 50, ylab = "# Observations", xlab = "Year", main = "Years")
  hist(as.numeric(data$MON), ylab = "# Observations", xlab = "months", main = "months (standardised to North Hemis)")
  #if(sum(!is.na(as.numeric(data$TIMEloc)) > 0)){
  #hist(as.numeric(data$TIMEloc), ylab = "# Observations", xlab = "Local Time of Day", main = "Time of Observations")
  #}
  hist(data$LATITUDE, ylab = "# Observations", xlab = "lat", main = "lat")

  
  print("Working on second set of histograms")
  par(mfrow = c(2,2))
  hist(log10(data$CHLO), ylab = "# Observations", xlab = expression(paste("log"[10], "(Chlo-a, mg m"^-3, ")")), 
       main = "Chlorophyll")
  hist(data$SST, ylab = "# Observations", xlab = "SST", main = "SST")    
  hist(data$BATHY, ylab = "# Observations", xlab = "Bathymetry", main = "Bathymetry")
  hist(log10(data$TOT_ABUND+min_val), ylab = "# Observations", xlab = expression(paste("log"[10], "(Abundance", ")")), 
       main = "# m^-3")
  
  print("Working on first set of abundance plots")
  par(mfrow = c(2,2))
  plot(data$BATHY, log10(data$TOT_ABUND+min_val), xlab = "Depth (metres)", 
       ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
       main = "Bathymetry", cex.lab = 1.3)
  plot(data$SST, log10(data$TOT_ABUND+min_val), xlab = "SST", 
       ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
       main = "SST", cex.lab = 1.3)
  plot(log10(data$CHLO), log10(data$TOT_ABUND+min_val), 
       xlab =expression(paste("log"[10], "(Chlo-a, mg m"^-3, ")")) , 
       ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
       main = "Chlorophyll", cex.lab = 1.3)
  plot(data$MESH, log10(data$TOT_ABUND+min_val), xlab = "Mesh Size (um)", 
       ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
       main = "Mesh Size", cex.lab = 1.3)
  
  print("Working on second set of abundance plots")
  par(mfrow = c(3,2))
  plot(data$MON, log10(data$TOT_ABUND+min_val), xlab = "month", 
       ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
       main = "month (North. Hemis.)", cex.lab = 1.3)
  if(sum(!is.na(as.numeric(data$TIMEloc))) > 0){
    plot(as.numeric(data$TIMEloc), log10(data$TOT_ABUND+min_val), 
         xlab = "Local Hour" , 
         ylab = expression(paste("log"[10], "(# m"^-3, ")")), 
         main = "Time of Day", cex.lab = 1.3)
  }
  plot(data$MID_Z, log10(data$TOT_ABUND+min_val),
       xlab = "Max. Tow Depth",
       ylab = expression(paste("log"[10], "(# m"^-3, ")")),
       main = "Sample Depth")
  plot(data$SST, log10(data$CHLO),xlab = "SST", 
       ylab = expression(paste("log"[10], "(Chlo-a, mg m"^-3, ")")) , 
       main = "CHLO & SST", cex.lab = 1.3)
  plot(data$MON, data$SST, 
       xlab = "month (North. Hemis.)" , 
       ylab = "SST", 
       main = "SST and month", cex.lab = 1.3)
  
  # Plot location of data for data
  print(paste("Working on global map of", zoo_name, "observations"))
  library(ggplot2)
  library(maptools)
  data(wrld_simpl)
  st_plot <- ggplot(data = data, aes(x = LONGITDE, y = LATITUDE, colour = 'black'))  + 
    geom_polygon(data=wrld_simpl, mapping=aes(x=long, y=lat, group = group), fill='grey', colour = "grey")+
    geom_point(colour = "black", size = 0.1)+
    theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
                       axis.title.y = element_blank(),
                       axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
                       panel.border = element_rect(colour = "black"),
                       panel.grid = element_blank(),
                       plot.margin = unit(c(0,0.2,0.,0.2), "cm"),
                       legend.position = "none",
                       plot.subtitle=element_text(size=12, color="black")) + coord_fixed(ratio = 1.)+ 
    xlab("") + ylab("")  +
    scale_x_continuous(limits = c(-180, 180), expand = c(0,0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0,0)) +
    ggtitle(label = paste(zoo_name, " observations", " (n = ", no_obs, ")", sep = ""))
  
  print(st_plot)
}
}
