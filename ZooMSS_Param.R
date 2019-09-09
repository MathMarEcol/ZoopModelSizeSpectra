###################################################
##### SET PARAMETERS
###################################################
rm()
params = list() # General and fish parameters

## Numerical integration stuff
params$dx = 0.1
params$tstepdays = 1
params$tmaxyears = 10


## log10(size ranges)
params$xmin = -12 # min phytoplankton size
params$xfmin = -1 # min fish size
params$xmax =  6 # max fish size
params$xs= 4 # age senescence mortality kicks in
params$xzs = -2
params$dx = 0.1
params$xpmax = -5
params$xzmin = -5 # max phytoplankton / min zooplankton size
params$xzmax = 0 # max zooplankton size
params$x = seq(params$xmin, params$xmax, params$dx)
params$y = seq(params$xmin, params$xmax, params$dx) # log equal size bins
params$xzminref = which(params$x==(params$xzmin)) # class at which dynamic equations begin
params$xzmaxref = which(params$x==params$xzmax)
params$xfminref = which(params$x==params$xfmin)
params$xmaxref = which(params$x==params$xmax)
params$w = 10^params$x
params$N = as.integer(365*params$tmaxyears/params$tstepdays)
params$dt = params$tmaxyears/params$N


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


params$D.0 <- 1 # theoretical minimum size (um)
params$m <- 0
params$sigma.z = 2.3
params$gamma.z = 875
params$alpha.z = 1.01
params$K.z = 0.7


