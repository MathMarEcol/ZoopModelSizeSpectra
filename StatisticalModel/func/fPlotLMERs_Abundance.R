fPlotLMERs_Abundance <- function (m1, Name, Transform) {
  # NOTE: To get CIs for lmer, need to use type = contrast (random effects cancel)
  Terms <- as.character(m1@call)[2] # Terms from the model so we can print blank if n.s.
  
  x11(width = 15, height = 6, title = Name)
  par(mfrow = c(2,5), mar = c(4,4,2,2))
  
  # if(grepl("Mesh", Terms, fixed = TRUE)) {
  #   visreg(m1, "Mesh", type = "contrast", scale = "response", xlab = "Mesh  (microns)")} else {
  #     plot.new()
  #   }
  # 
  # if(grepl("Tow", Terms, fixed = TRUE)) {
  #   visreg(m1, "Tow", rug = FALSE, type = "contrast", scale = "response", xlab = "Tow")} else {
  #     plot.new()
  #   }
  
  if(grepl("Chl", Terms, fixed = TRUE)) {
    visreg(m1, "Chl", type = "contrast", scale = "response", xlab = "Chl-a (mg/m3)")} else {
      plot.new()
    }
  if(grepl("Bathy", Terms, fixed = TRUE)) {
    visreg(m1, "Bathy", type = "contrast", scale = "response", xlab = "Bathy (m)")} else {
      plot.new()
    }
  
  if(grepl("SST", Terms, fixed = TRUE)) {
    visreg(m1, "SST", type = "contrast", scale = "response", xlab = "SST (ºC)")} else {
      plot.new()
    }
  
  if(grepl("HarmDOY", Terms, fixed = TRUE)) {
    visreg(m1, "HarmDOY", type = "contrast", scale = "response", xlab = "Day of year")} else {
      plot.new()
    }
  
  visreg2d(m1, yvar = "HarmDOY", xvar = "SST", scale = "response", plot.type = "persp", theta = 45, phi = 10, r = 100, 
           ticktype = "detailed", xlab = "\nSST (ºC)", ylab = "\nDay of Year", zlab = "\nlog10(Abundance)", color = "deepskyblue2")
  
  if(grepl("Mid_Z", Terms, fixed = TRUE)) {
      visreg(m1, "Mid_Z", type = "contrast", scale = "response", xlab = "Depth")} else {
      plot.new()
    }
  
  if(grepl("HarmHour", Terms, fixed = TRUE)) {
    visreg(m1, "HarmHour", type = "contrast", scale = "response", xlab = "Time of Day")} else{
      plot.new()
    }

  if(grepl("exp(-Mid_Z/1000) * fHarmonic(HarmHour, k = 1)", Terms, fixed = TRUE) | grepl("exp(-Mid_Z/1000):fHarmonic(HarmHour, k = 1)", Terms, fixed = TRUE)) {
    visreg(m1, "HarmHour", by = "Mid_Z", breaks = c(0, 50, 100, 300, 500, 1000), xlab = "Time of day", ylab = "log10(Abundance)", 
           type = "contrast", scale = "response", overlay = TRUE, rug = 0, strip.names = c("Depth=0","Depth=50", "Depth=100", "Depth=300", "Depth=500", "Depth=1000"))} else {
             plot.new()
           }
  
  # dev.copy2eps(file = paste0("Figures/", Name, ".eps"))
  # Shows the surface in colour, but does not have CIs on Depth*TOD
  dev.print(pdf, paste0("Figures/", Name, "_", Transform, ".pdf"))
}