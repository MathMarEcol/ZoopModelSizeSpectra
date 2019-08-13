fSummariseLMERs <- function (m1, Name) {
  # Summary, r2 and diagnostic plots
  r2 <- r.squaredGLMM(m1)
  
  # # Diagnostic plots
  # x11(width = 8, height = 5, title = paste0(Name, "Diagnostics"))
  # par(mfrow = c(1,3))
  # qqnorm(residuals(m1))
  # qqline(residuals(m1), col = "blue")
  # plot(ranef(m1)$Longhurst, main = "Longhurst")
  
  # if(Name == "Euphausiids")
  #   plot(ranef(m1)$Project, main = "Project") else 
  #     plot(ranef(m1)$Transect, main = "Transect")
  
  # Diagnostic plots
  x11(width = 8, height = 5, title = paste0(Name, "Diagnostics"))
  par(mfrow = c(1,3))
  qqnorm(residuals(m1))
  qqline(residuals(m1), col = "blue")
  plot(ranef(m1), main = "Mesh_Grid_Tow")
  
  return(r2)
}