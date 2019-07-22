# https://stats.stackexchange.com/questions/131106/predicting-with-random-effects-in-mgcv-gam

library(tidyverse)
library(raster)
library(rasterImage)
library(sf)
library(ggpubr)
library(splines)

source("fHarmonic.R")
source("fPlotGlobalGAMM.R")

####################################################
NthMonth <- c("01January", "02February","03March","04April","05May","06June",
              "07July","08August","09September","10October","11November","12December")
              
# Select the Taxa you want to plot
taxa <- c("OmniCopepods", "CarnCopepods", "Chaetognaths",
          "Salps", "Euphausiids", "Larvaceans", "Jellyfish")

taxa <- c("Euphausiids", "Jellyfish")

df <- data.frame(HarmHour = 1, Mesh = 100, Tow = "S1",
           Mid_Z = 10, Gear = 191, Longhurst = "ANTA", Type = "CPR", 
           Project = "IMOS-CPR", Transect = "AA")

# Set the max values to be removed in rasterbrick
maxm <- tibble(Bathy = 6000, SST = 30, Chl = 10)

##############################################################################

# Predict abundance using the predict function in raster and plot
tframe <- "annual"

for (i in 1:length(taxa)){
  if (taxa[i]=="Euphausiids"){
    excl <- c("s(Longhurst)", "s(Project)")} else{ # Zero out Project for Euphausiids
    excl <- c("s(Longhurst)", "s(Transect)")} # Otherwise Transect for everything else
  myplots <- fPlotGlobalGAMM(NthMonth,taxa[i],df,excl,tframe, maxm)
}
