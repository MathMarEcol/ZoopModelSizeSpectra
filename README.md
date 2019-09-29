# Model code for Zooplankton Model of Size Spectrum (ZooMSS v1)

This repo holds the code for Zooplankton Model of Size Spectrum (ZooMSS v1) as published in Heneghan et al. 2016 (pdf attached).   

The size-spectrum model consists of a size-spectrum comprised of three communities: phytoplankton, zooplankton, and fish. The phytoplankton component covers the smallest size classes and is held constant as a background resource spectrum for zooplankton. Size-dependent processes of growth and mortality drive the zooplankton and fish components.  

Heneghan RF, Everett JD, Blanchard JL and Richardson AJ (2016) Zooplankton Are Not Fish: Improving Zooplankton Realism in Size-Spectrum Models Mediates Energy Transfer in Food Webs. *Front. Mar. Sci.* **3**:201. doi: 10.3389/fmars.2016.00201
https://doi.org/10.3389/fmars.2016.00201  

## SCRIPTS 
ZooMSS_Shell.R: You run the model from this script.    
ZooMSS_Parameters.R: Parameters for the model.   
ZooMSS_Stability: Code for stability analysis via Newton-Raphson method.    
ZooMSS_PDE: The model code.    
