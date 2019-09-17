## Shell to run ZooMSS v1 (see Heneghan et al. 2016: https://doi.org/10.3389/fmars.2016.00201)
## Author: Ryan Heneghan
## Last updated: September 2019

rm(list = ls())

source("ZooMSS_Param.R") # Parameter list
source("ZooMSS_PDE.R") # PDE implementation
source("ZooMSS_Stability.R") # Stability analysis code

attach(params)
## Set up initial value vectors for zooplankton and fish
init.spec = a*w^b
N.p = init.spec*I(x<x[xzminref]) # Phytoplankton spectrum (held constant at the moment)
N.z.init = init.spec[x>=x[xzminref] & x<=x[xzmaxref]] # Zooplankton spectrum
N.f.init = init.spec[x>=x[xfminref]] # Fish spectrum
init.dist = c(N.z.init, N.f.init)
detach(params)

## Run ZooMSS model, test = 1 if you want a progress bar and live plot
out <- ZooMSS_PDE(state = init.dist, parms = params, test = 1)

## Plot average size spectra from last 50% of model runs
slope.matrix0 <- matrix(unlist(out[1]), nrow = length(init.dist), ncol = params$N) # Extract abundances
ave.size.spectra0 = rowMeans(slope.matrix0[,c((0.5*params$N):(params$N-1))]) # Calcualte average spectra from final 50%

# Create plot
plot(log10(params$w), log10(N.p), col = 'green', lwd = 2, type = 'l', ylim = c(log10(min(ave.size.spectra0)),log10(N.p[1])), 
     ylab = expression(paste('log'[10],'(Abundance # m'^{-3},')')), xlab = expression(paste('log'[10], '(Body Weight, g)')))
lines(params$x[params$x>=params$x[params$xzminref] & params$x<=params$x[params$xzmaxref]], log10(ave.size.spectra0[1:51]), lwd = 2, col = 'red')
lines(params$x[params$x>=params$x[params$xfminref]], log10(ave.size.spectra0[52:122]), lwd = 2, col = 'blue')
legend('bottomleft', legend = c('Phytoplankton', 'Zooplankton', 'Fish'), fill = c('green', 'red', 'blue'), bty = 'n')
