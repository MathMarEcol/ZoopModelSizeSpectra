###################################################
##### SET PARAMETERS FOR ZooMSS v1
###################################################
rm()
params = list() # General and fish parameters

## Numerical integration stuff
params$dx = 0.1 # width of log size bings
params$tstepdays = 1 # time step in days
params$tmaxyears = 10 # length of simulation

## log10(size ranges)
params$xmin = -12 # min phytoplankton size
params$xfmin = -1 # min fish size
params$xmax =  6 # max fish size
params$xs= 4 # size at which senescence mortality for fish kicks in
params$xzs = -2 # size at which senescence mortality for zooplankton kicks in
params$xpmax = -5 #  max phytoplankton size
params$xzmin = -5 # min zooplankton size
params$xzmax = 0 # max zooplankton size
params$x = seq(params$xmin, params$xmax, params$dx) # log equal size bins
params$y = seq(params$xmin, params$xmax, params$dx) # log equal size bins
params$xzminref = which(params$x==(params$xzmin)) # size bin at which dynamic equations begin
params$xzmaxref = which(params$x==params$xzmax) # size bin at which zooplankton spectrum ends
params$xfminref = which(params$x==params$xfmin) # size bin at which fish spectrum beings 
params$xmaxref = which(params$x==params$xmax) # max size bin
params$w = 10^params$x # size classes in normal space
params$N = as.integer(365*params$tmaxyears/params$tstepdays) # number of time steps
params$dt = params$tmaxyears/params$N # time step size

## Parameter values
params$sigma.f = 2.3 
#params$r = 0.25
params$alpha.f = 0.82
params$gamma.f = 640
params$K.f = 0.6
params$rho = 0.1
params$S.0 = 0.04
params$s = -0.25
params$a = 0.017
params$b = -1
params$k.sm = 0.2
params$k.zsm = 0.05
params$p.s = 0.3
params$p.zs = 1.2
params$beta.fish = 100
params$beta.zoo = 100
params$r.0 = 1
params$p = 0.75
params$kappa.p = 5e-2
params$lambda = 1
params$r = 1
params$n = 0.75

## Zooplankton specific parameters
params$D.0 = 1 
params$m = 0
params$sigma.z = 2.3
params$gamma.z = 875
params$alpha.z = 1.01
params$K.z = 0.7


