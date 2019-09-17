## Numerical implementation of ZooMSS v1 (see Heneghan et al. 2016: https://doi.org/10.3389/fmars.2016.00201)
## Author: Ryan Heneghan
## Last updated: September 2019

ZooMSS_PDE <- function(state, parms, test){
  with(as.list(c(state, parms)),{
    
    ################### CONSTANT STUFF ################### 
    ## Beta values for zooplankton and fish
    beta = matrix(0,2,length(x))
    beta[1,c(xfminref:xmaxref)] = log(beta.fish) # beta for fish
    
    # Functions to calculate zooplankton beta
    D.z <- 2*(3*w[c((xzminref-1):xzmaxref)]*1e12/(4*pi))^(1/3)# convert body mass g to ESD (um)
    
    beta.z <- function(m.){
      betaz =  log((exp(0.02*log(D.z/D.0)^2 - m. + 1.832))^3)
      return(betaz)
    }
    beta[2,c((xzminref-1):xzmaxref)] <-  beta.z(m)
    
    q1 = matrix(NA, length(x), length(y))
    lx = ly = log(w)
    for (i in 1:length(y)) { q1[,i] = ly[i] - lx}
    
    betafish.mat <- matrix(beta[1,], nrow = length(x), ncol = length(y) , byrow = TRUE)
    betazoo.mat <- matrix(beta[2,], nrow = length(x), ncol = length(y), byrow = TRUE) 
    
    phi.f = function(qmat,beta.mat,sigma, A, alpha){ # feeding kernel function
      qtemp = qmat - beta.mat
      phi=ifelse(beta.mat != 0,exp(-(qtemp)^2/(2*sigma*sigma))/(sigma*sqrt(2*pi)),0)
      gphi = phi # growth kernel matrix
      
      ## Multiply by search rate
      search = matrix(A*10^(alpha*x), nrow = length(x), ncol = length(y), byrow = TRUE)
      mphi = phi*search # mortality kernel matrix (prey are rows, predators are columns)
      gphi = t(mphi) # growth kernel matrix (predators are rows, prey are columns)
      return(list(gphi,mphi))
    }
    
    ## Weight diff matrices
    
    wgdiff <- (10^-x)%*%t(10^x)/log(10)
    wddiff <- (10^-x)%*%t(10^(2*x))/log(10)
    
    ## Growth efficiency rates
    effic.fish = matrix(I(x>=x[xfminref])*(1-(w/w[xmaxref])^(r-n)), nrow = length(x), ncol = length(y))
    
    effic.zoo = matrix(I(x >= x[xzminref-1] & x <= x[xzmaxref])*(1-(w/w[xzmaxref])^(r-n)), nrow = length(x), ncol = length(y))
    effic.zoo[is.nan(effic.zoo)] = 0
    
    ## Simpson's Rule vector for integration
    simp <- array(1, dim = length(x))
    simp[c(seq(2,length(x)-1,2))] = 4
    simp[c(seq(3,length(x)-2,2))] = 2
    
    sm <- matrix(simp, nrow = length(x), ncol = length(x), byrow = TRUE)*(dx/3)
    
    # Calculate matrices for all constant parts of growth, diffusion and death integrals
    gphi.f = K.f*matrix(unlist(phi.f(q1, betafish.mat, sigma.f, gamma.f, alpha.f)[1]), nrow = length(x), ncol = length(y))*wgdiff*sm # fish growth integral matrix
    dphi.f = (K.f^2/2)*matrix(unlist(phi.f(q1, betafish.mat, sigma.f, gamma.f, alpha.f)[1]), nrow = length(x), ncol = length(y))*wddiff*sm # fish diffusion integral matrix
    mphi.f = matrix(unlist(phi.f(q1, betafish.mat, sigma.f, gamma.f, alpha.f)[2]), nrow = length(x), ncol = length(y))*sm # fish mortality integral matrix
    
    gphi.z = K.z*matrix(unlist(phi.f(q1, betazoo.mat, 1, gamma.z, alpha.z)[1]), nrow = length(x), ncol = length(y))*wgdiff*sm # fish growth integral matrix
    dphi.z = (K.z^2/2)*matrix(unlist(phi.f(q1, betazoo.mat, sigma.z, gamma.z, alpha.z)[1]), nrow = length(x), ncol = length(y))*wddiff*sm # fish diffusion integral matrix
    mphi.z = matrix(unlist(phi.f(q1, betazoo.mat, sigma.z, gamma.z, alpha.z)[2]), nrow = length(x), ncol = length(y))*sm # fish mortality integral matrix
    
    ## Time steps vector for ode
    time.steps <- seq(0, tmaxyears, dt)
    
    ####################################################	
    
    ################ DYNAMIC EQUATIONS #################
    N.z.st = c(array(0,c(1,(xzminref-1))), state[1:(xzmaxref-xzminref + 1)], array(0, c(1, c(xmaxref-xzmaxref))))
    N.f.st = c(array(0,c(1,xfminref-1)), state[(xzmaxref-xzminref + 2): length(state)])
    
    ## Storage arrays
    # Matrices for recording phytoplankton, zooplankton, fish and community size spectra
    N.z = N.f = N.c = array(0, c(length(x), N)) 
    state.new = array(0, c(length(state), N))
    
    # Matrices for keeping track of ingested food
    F.z = F.f = G.z = G.f = array(0, c(length(x), N))
    
    # Matrices for keeping track of mortality
    M.z = M.f = array(0, c(length(x), N)) 
    
    N.z[,1] = N.z.st
    N.f[,1] = N.f.st
    N.c[,1] = init.spec
    
    if(test == 1){
      pb <- txtProgressBar(min = 0, max = N, style = 3)
      plot(x, log10(N.c[,1]), type="l", xlim=c(xmin,xmax))
    }
    
    state.new[,1] = state
    
    for(i in 1:(N-1)){
      
      if(test == 1){setTxtProgressBar(pb,i)}
      
      #### Calculate growth, mortality and diffusion integrals
      N.pmat = matrix(N.p, nrow = length(x), ncol = length(x), byrow = TRUE)
      N.zmat = matrix(N.z[,i], nrow = length(x), ncol = length(x), byrow = TRUE)
      N.fmat = matrix(N.f[,i], nrow = length(x), ncol = length(x), byrow = TRUE)
      
      # Background and senescence mortality
      BM.z = (S.0*w^(s) + k.zsm*10^(p.zs*(x-xzs)))*I(x>=x[xzminref] & x<=x[xzmaxref])
      BM.f = (S.0*w^(s) + k.sm*10^(p.zs*(x-xs)))*I(x>=x[xfminref]) 
      
      ## Zooplankton integrals
      growth.z = (rowSums((gphi.z*N.pmat)) + rowSums((gphi.z*N.zmat)) + rowSums((gphi.z*N.fmat)))
      mort.z = rowSums((mphi.z*N.pmat)) + rowSums((mphi.z*N.zmat)) + rowSums((mphi.z*N.fmat)) + BM.z
      diff.z = (rowSums((dphi.z*N.pmat)) + rowSums((dphi.z*N.zmat)) + rowSums((dphi.z*N.fmat)))
      
      ## Fish integrals
      growth.f = (rowSums((gphi.f*N.pmat)) + rowSums((gphi.f*N.zmat)) + rowSums((gphi.f*N.fmat)))
      mort.f = rowSums((mphi.f*N.pmat)) + rowSums((mphi.f*N.zmat)) + rowSums((mphi.f*N.fmat)) + BM.f
      diff.f = (rowSums((dphi.f*N.pmat)) + rowSums((dphi.f*N.zmat)) + rowSums((dphi.f*N.fmat)))
      
      ### Store growth and mortality
      F.z[,i] = growth.z/K.z
      F.f[,i] = growth.f/K.f
      G.z[,i] = growth.z
      G.f[,i] = growth.f
      M.z[,i] = mort.z
      M.f[,i] = mort.f
      
      ########################################################
      
      ####### ZOOPLANKTON
      #### Solve Mvf first
      N.z.iter = array(0, c(1, (xzmaxref - xzminref + 1)))
      G.0 = G.z[xzminref-1,i]
      A.z.iter = 1 + dt/dx*G.z[(xzminref:xzmaxref),i] + dt*M.z[(xzminref:xzmaxref),i]
      B.z.iter = c(dt/dx*G.0, dt/dx*G.z[(xzminref:(xzmaxref-1)),i])
      
      N.z.iter = (N.z[c(xzminref:xzmaxref),i] + c(N.p[xzminref-1], N.z[xzminref:(xzmaxref-1),i])*B.z.iter)/A.z.iter
      
      ### Solve MvF with diffusion
      A.z = c(dt/dx*G.0, dt/dx*G.z[(xzminref:(xzmaxref-1)),i]) +
        c(0, diff.z[(xzminref:(xzmaxref-2))]*dt/(2*dx^2), 0)
      B.z = 1 + dt/dx*G.z[(xzminref:xzmaxref),i] + dt*M.z[(xzminref:xzmaxref),i] + 
        c(0, dt/(dx^2)*diff.z[(xzminref +1):(xzmaxref - 1)] ,0)
      C.z = c(0, dt/(2*dx^2)*diff.z[(xzminref+2): xzmaxref]*N.z.iter[3:length(N.z.iter)], 0)
      N.z[(xzminref:xzmaxref),i+1] = (N.z[c(xzminref:xzmaxref),i] + 
                                        c(N.p[xzminref-1], N.z[xzminref:(xzmaxref-1),i])*A.z + C.z)/B.z
      
      ####### FISH
      #### Solve Mvf first
      N.f.iter = array(0, c(1, (xmaxref - xfminref + 1)))
      G.0 = G.z[xfminref-1,i]
      A.f.iter = 1 + dt/dx*G.f[(xfminref:xmaxref),i] + dt*M.f[(xfminref:xmaxref),i]
      B.f.iter = c(dt/dx*G.0, dt/dx*G.f[(xfminref:(xmaxref-1)),i])
      
      N.f.iter = (N.f[c(xfminref:xmaxref),i] + c(N.z[xfminref-1], N.f[xfminref:(xmaxref-1),i])*B.f.iter)/A.f.iter
      
      ### Solve MvF with diffusion
      A.f = c(dt/dx*G.0, dt/dx*G.f[(xfminref:(xmaxref-1)),i]) +
        c(0, diff.f[(xfminref:(xmaxref-2))]*dt/(2*dx^2), 0)
      B.f = 1 + dt/dx*G.f[(xfminref:xmaxref),i] + dt*M.f[(xfminref:xmaxref),i] + 
        c(0, dt/(dx^2)*diff.f[(xfminref +1):(xmaxref - 1)] ,0)
      C.f = c(0, dt/(2*dx^2)*diff.f[(xfminref+2): xmaxref]*N.f.iter[3:length(N.f.iter)], 0)
      
      N.f[(xfminref:xmaxref),i+1] = (N.f[c(xfminref:xmaxref),i] + 
                                       c(N.z[xfminref-1], N.f[xfminref:(xmaxref-1),i])*A.f + C.f)/B.f
      N.f[xfminref, i+1] = N.z[xfminref, i+1]
      
      ####### COMMUNITY
      state.new[,i+1] = c(N.z[(xzminref:xzmaxref),i+1], N.f[(xfminref:xmaxref),i+1])
      
      ####### LIVE PLOT
      if(test == 1){
        N.c[1:(xzminref-1),i+1] = N.p[1:(xzminref-1)]
        N.c[xzminref:(xfminref-1), i+1] = N.z[xzminref:(xfminref-1),i+1] 
        N.c[xfminref:xmaxref, i+1] = N.f[xfminref:xmaxref,i+1]
        
        r = rainbow(N, s=1, v=1, start=0, end=max(1,N - 1)/N)
        lines(x, log10(N.z[,i]), type="l", col=r[i], cex=1.2)
        lines(x, log10(N.f[,i]), type="l", col =r[i], cex = 1.2)
      }
      
    } # End MvF-D loop
    return(list(state.new, F.z, F.f, M.z, M.f)) # state.new is saved abundances of zoo and fish
                                                # F.z is feeding rate of zooplankton
                                                # F.f is feeding rate of fish
                                                # M.z is mortality rate of zooplankton
                                                # M.f is mortality rate of fish
    
  }) # End with(as.list())...
  
} # End ZooMSS_PDE