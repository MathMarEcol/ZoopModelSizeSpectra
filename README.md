# Model code for Zooplankton Model of Size Spectrum (ZooMSS)


## Overview of ZooMSS

The Zooplankton Model of Size Spectra (ZooMSS) is a functional size-spectrum model of the marine ecosystem (following Heneghan et al. 2016, ZooMSSv1) to resolve phytoplankton, nine zooplankton functional groups (heterotrophic flagellates and ciliates, omnivorous and carnivorous copepods, larvaceans, euphausiids, salps, chaetognaths and jellyfish) and three size-based fish groups. Zooplankton functional groups are resolved using their body-size ranges, size- based feeding characteristics and carbon content, and the zooplankton community emerges from the model across global environmental gradients, depending on the functional traits of the different groups. 

We developed the Zooplankton Model of Size Spectra version 2 (ZooMSSv2) based on the prototype of Heneghan et al. (2016). ZooMSSv2 uses the functional size-spectrum framework (Blanchard et al., 2017) to resolve the body size ranges, size-based feeding characteristics and carbon content of nine zooplankton groups and three fish groups. The model is run on 5° grid squares across the global ocean. For each region, the model is forced with the long- term mean satellite sea surface temperature and chlorophyll a from MODIS-Aqua.

ZooMSSv2 represents the marine ecosystem as three communities: phytoplankton, zooplankton and fish. The zooplankton community consists of nine of the most abundant zooplankton groups, and the fish community was made up of a small, medium and large group. Dynamics of the phytoplankton are not explicitly resolved in the model, rather the mean size structure of the phytoplankton community is estimated directly from satellite chlorophyll a observations (Brewin et al., 2010; Barnes et al., 2011; Hirata et al., 2011). Abundances of the zooplankton and fish communities are driven by size-dependent processes of growth and mortality, with the temporal dynamics of each functional group governed by separate second- order McKendrick-von Foerster equations.


## How to run ZooMSS using `R`

### Model Structure

The model is split into 5 main files.

1. Setup_RUN_NAME.R - Set up your own specific simulation. This is the only file you usually need to edit.

2. fZooMSS_Model.R - The top-level model function which seuqentially runs the following 3 functions, and saves the output.

    + fZooMSS_Params.R - Set all the parameter values for the model.
  
    + fZooMSS_Setup.R - Do all calculations which can be done in advance of the time loop.
  
    + fZooMSS_Run.R - Contains the time-dependent loop.

### Input files

ZooMSS requires two input files:   

1. TestGroups.csv which contains all the taxa-specific parameter values, including size ranges.   

2. And RDS file which contains all environmental data for each grid cell. This includes chlorophyll a (chlo), sea surface temperature (sst), bathymetry (bathy), phytoplankton slope and intercept, and latitude/longitude which is useful for plotting.   


### Running
To run the model, you simply need to provide the relevent information in Setup_RUN_NAME.R and run the script.

### Plotting
Besides writing your own plotting routines, you can now easily plot the overall species-resolved size-spectra and the PPMRs.
Run using:    
* `gg <- ZooMSS_Plot_SizeSpectra(dat)`    
* `gg <- ZooMSS_Plot_PPMR(dat)`    

If you have saved the timesteps you can plot the timeseries:    
* `gg <- ZooMSS_Plot_AbundTimeSeries(dat)`  
* `gg <- ZooMSS_Plot_GrowthTimeSeries(dat)`   
* `gg <- ZooMSS_Plot_PredTimeSeries(dat)` 
* `gg_list <- fZooMSS_Plot_RuntimeAveraging(file_name)`   
* `fZooMSS_Plot_Wavelet(file_name, species_name)`   

where `dat` is the model output list, and `gg` is a ggplot output from the plotting routines.

### ZooMSS Dashboard
Assuming you have saved all the time-steps (using `SaveTimeSteps <- TRUE`) you can explore the output using our newly created **ZooMSS Dashboard**. To access, go here: https://jaseeverett.shinyapps.io/ZooMSS_Dashboard/. Please drop us a line with any comments or ideas.

## Publications

1. Heneghan, R.F., Everett, J.D., Blanchard, J.L., Richardson, A.J., 2016. Zooplankton Are Not Fish: Improving Zooplankton Realism in Size-Spectrum Models Mediates Energy Transfer in Food Webs. Front. Mar. Sci. 3, 1–15. https://doi.org/10.3389/fmars.2016.00201

2. Heneghan, R.F., Everett, J.D., Sykes, P., Batten, S.D., Edwards, M., Takahashi, K., Suthers, I.M.,  Blanchard, J.L., Richardson, A.J., 2020. A global size-spectrum model of the marine ecosystem that resolves zooplankton composition. Ecological Modelling 435, 109265. https://doi.org/10.1016/j.ecolmodel.2020.109265


## Getting Help

If you are running the model come across a problem, raise an issue on GitHub: https://github.com/MathMarEcol/ZoopSizeSpectraModel/issues

If you find some errors, or just want to improve the model, we'd love you to go ahead and make the changes and then submit a pull request for us to check and approve.


