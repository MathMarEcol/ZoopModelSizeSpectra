## Newton Raphson stability analysis code for ZooMSS v1 (see Heneghan et al. 2016: https://doi.org/10.3389/fmars.2016.00201)
## Author: Ryan Heneghan
## Last updated: September 2019

########## NEWTON RAPHSON
## du.calc:
## u_initial is the starting size spectrum

du.calc <- function(u_initial){
# Run the PZF model, output is a matrix of size spectrum for each time step
# the model is run (usually just N = 2, so only 2 size spectra)
out <- ZooMSS_PDE(state = u_initial, parms = params, test = 0)[1]

# Take the last time step's size spectrum as the new size spectrum to calculate
# du for each size class
oldvals <- u_initial
newvals <- matrix(unlist(out), nrow = length(oldvals), ncol = params$N)[,params$N]

# Calculate du
du <- as.vector((newvals - oldvals)/params$dt)
return(du)
}

## Calculate Jacobian Matrix
jac.calc <- function(u.init, du){
	jac = matrix(0, length(du), length(du))
	
	for(j in 1:length(du)){
	u_temp = u.init  
	temp=u.init[j] #Select size class
    
  h=(1e-4)*abs(temp)	# Calculate perturbation for current size class, depends on
                      # abundance in current size class being evaluated (temp)
	
	if(h == 0.0){h = 1e-4} # If there's nothing in the current size class, fix the
                         # perturbation to a minimum (1e-4)

    u_temp[j]=(temp+h) #Increment variable

 	  h=u_temp[j]-temp  #Error rounding check step

    du.new <- du.calc(u_temp) #Calculate new rhs values

    u_temp[j]=temp #Restore variable value

    jac[,j]=((du.new-du)/h) #Calculate derivatives for Jacobian Matrix
  }
  return(jac)
  }



Newton <- function(u_initial, ntrials = 20, tol = 1e-4){
	u_temp<-u_initial
  erf = 1 # Starting value for convergence function
  i = 0
  while(erf > tol && i < ntrials){ # This runs until convergence
                                   # or we run out of trials
    
    du <- du.calc(u_temp) ## Calculate current du's

    jac.c <- jac.calc(u_temp, du) ## Calculate current Jacobian

    u_temp <- as.vector(u_temp + solve(jac.c,-du, tol=1e-25)) # Update u_temp
    
   
    erf<-sum(abs(du))  # Check sum of absolute values of du 
                       # (if erf gets smaller, the size spectrum is settling down)
    
    i = i + 1
  }
  return(max(Re(eigen(jac.c)$values))) ## Return lambda_max
}
