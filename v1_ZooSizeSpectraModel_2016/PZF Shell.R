
setwd("/Users/ryanheneghan/Library/Mobile Documents/com~apple~CloudDocs/PhD/Writing/Frontiers Manuscript/R Code")

rm(list = ls())
source("PZF_Parameters.R")
source("ODE Algorithm.R")
source("PZF Stability.R")
################## 

attach(params)
## Set up initial value vectors for zooplankton and fish
init.spec = a*w^b
N.p = init.spec*I(x<x[xzminref]) # Phytoplankton spectrum (held constant at the moment)
N.z.init = init.spec[x>=x[xzminref] & x<=x[xzmaxref]] # Zooplankton spectrum
N.f.init = init.spec[x>=x[xfminref]] # Fish spectrum
init.dist = c(N.z.init, N.f.init)
detach(params)

params$N <- 365*20

########## VALUE RANGES #############
sigma.slide <- seq(0.9, 2.5, 0.1)
m.slide <- seq(-3, 2, 0.1)
K.slide <- seq(0.3,0.9,0.05)
l.m <- length(m.slide)
l.s <- length(sigma.slide)
l.K <- length(K.slide)

######### STORAGE MATRICES ###########
store.names <- c("Biomass", "P:B", "TP", "CV", "EE")
sigma.store <- matrix(0, nrow = l.s, ncol = 5)
m.store <- matrix(0, nrow = l.m, ncol = 5)
K.store <- matrix(0, nrow = l.K, ncol = 5)
colnames(sigma.store) <- store.names
rownames(sigma.store) <- sigma.slide
colnames(m.store) <- store.names
rownames(m.store) <- m.slide
colnames(K.store) <- store.names
rownames(K.store) <- K.slide

zoop.range = 1:(params$xzmaxref - params$xzminref + 1)
zoop.weights = params$w[params$xzminref:params$xzmaxref]
fish.weights = params$w[params$xfminref:params$xmaxref]	
	
N.runs = l.m + l.s + l.K


###### BASE MODEL
########################################################################################
# Base Model
params$K.z = params$K.f
params$sigma.z = params$sigma.f
params$gamma.z = params$gamma.f
params$alpha.z = params$alpha.f
params$beta.zoo = params$beta.fish

out <- PZF.solve(state = init.dist, parms = params, test = 1, zoo.beta = 0)


slope.matrix0 <- matrix(unlist(out[1]), nrow = length(init.dist), ncol = params$N)
ave.size.spectra0 = rowMeans(slope.matrix0[,c((0.5*params$N):(params$N-1))])
feeding.fish.matrix <- matrix(unlist(out[3]), nrow = length(params$x), ncol = params$N)
mort.fish.matrix <- matrix(unlist(out[5]), nrow = length(params$x), ncol = params$N)
ave.mort.fish = rowMeans(mort.fish.matrix[(params$xfminref:params$xmaxref), 
											c((0.5*params$N):(params$N-1))])
ave.thru.fish = rowMeans(feeding.fish.matrix[(params$xfminref:params$xmaxref),
												c((0.5*params$N):(params$N-1))])		
pbf0 = sum(ave.mort.fish*ave.size.spectra0[-zoop.range]*fish.weights*params$dx)/(fish.0)			
tpf0 = sum(ave.thru.fish*ave.size.spectra0[-zoop.range]*fish.weights*params$dx)										
zoo.0 = sum(ave.size.spectra0[zoop.range]*zoop.weights*params$dx)
									# Average biomass density of zooplankton
fish.0 = sum(ave.size.spectra0[-zoop.range]*fish.weights*params$dx)
									# Average biomass density of fish

ee0 = fish.0/zoo.0
								
######
params$K.z = 0.7
params$gamma.z = 875
params$alpha.z = 1.01

for(i in 1:N.runs){
	print(i/N.runs)
	
	if(i <= l.m){
	params$m <- m.slide[i]
	params$K.z <- 0.7
	params$sigma.z <- 1.75
		
	out <- PZF.solve(state = init.dist, parms = params, test = 1, zoo.beta = 1)
			
	slope.matrix <- matrix(unlist(out[1]), nrow = length(init.dist), ncol = params$N)
	feeding.zoo.matrix <- matrix(unlist(out[2]), nrow = length(params$x), ncol = params$N)
	feeding.fish.matrix <- matrix(unlist(out[3]), nrow = length(params$x), ncol = params$N)
	mort.zoo.matrix <- matrix(unlist(out[4]), nrow = length(params$x), ncol = params$N)
	mort.fish.matrix <- matrix(unlist(out[5]), nrow = length(params$x), ncol = params$N)	
	
				
	##### FISH BIOMASS
	slope.store = rowMeans(slope.matrix[,c((0.5*params$N):(params$N-1))])
	zoo.biom = sum(slope.store[zoop.range]*zoop.weights*params$dx)
	fish.biom = sum(slope.store[-zoop.range]*fish.weights*params$dx)
			
	##### P:B RATIO
	ave.mort.fish = rowMeans(mort.fish.matrix[(params$xfminref:params$xmaxref), 
																	c((0.5*params$N):(params$N-1))])
																	
	prod.store = sum(ave.mort.fish*slope.store[-zoop.range]*fish.weights*params$dx)

	pb.store = prod.store/fish.biom

	##### TOTAL THROUGHPUT
	ave.thru.fish = rowMeans(feeding.fish.matrix[(params$xfminref:params$xmaxref),																					 c((0.5*params$N):(params$N-1))])

	tp.store = sum(ave.thru.fish*slope.store[-zoop.range]*fish.weights*params$dx)
																		
			
	##### COEFFICIENT OF VARIATION
	cv = sd(colMeans(fish.weights*slope.matrix[-zoop.range, c((0.5*params$N):(params$N-1))]*params					$dx))/
				mean(colMeans(fish.weights*slope.matrix[-zoop.range, c((0.5*params$N):(params$N-1))]*params					$dx))						
			
	##### ECOLOGICAL EFFICIENCY
	ecol.effic = (fish.biom/zoo.biom)/ee0
	
	m.store[i,1] = fish.biom
	m.store[i,2] = pb.store
	m.store[i,3] = tp.store
	m.store[i,4] = cv
	m.store[i,5] = ecol.effic
	}
	
	if(i %in% ((l.m+1):(l.m+l.s))){
	params$sigma.z <- sigma.slide[(i-l.m)]
	params$m <- 0
	params$K.z <- 0.7
	i.mod <- (i-l.m)	
	
	out <- PZF.solve(state = init.dist, parms = params, test = 1, zoo.beta = 1)
			
	slope.matrix <- matrix(unlist(out[1]), nrow = length(init.dist), ncol = params$N)
	feeding.zoo.matrix <- matrix(unlist(out[2]), nrow = length(params$x), ncol = params$N)
	feeding.fish.matrix <- matrix(unlist(out[3]), nrow = length(params$x), ncol = params$N)
	mort.zoo.matrix <- matrix(unlist(out[4]), nrow = length(params$x), ncol = params$N)
	mort.fish.matrix <- matrix(unlist(out[5]), nrow = length(params$x), ncol = params$N)	
		##### FISH BIOMASS
	slope.store = rowMeans(slope.matrix[,c((0.5*params$N):(params$N-1))])
	zoo.biom = sum(slope.store[zoop.range]*zoop.weights*params$dx)
	fish.biom = sum(slope.store[-zoop.range]*fish.weights*params$dx)
			
	##### P:B RATIO
	ave.mort.fish = rowMeans(mort.fish.matrix[(params$xfminref:params$xmaxref), 
																	c((0.5*params$N):(params$N-1))])
																	
	prod.store = sum(ave.mort.fish*slope.store[-zoop.range]*fish.weights*params$dx)

	pb.store = prod.store/fish.biom

	##### TOTAL THROUGHPUT
	ave.thru.fish = rowMeans(feeding.fish.matrix[(params$xfminref:params$xmaxref),																					 c((0.5*params$N):(params$N-1))])

	tp.store = sum(ave.thru.fish*slope.store[-zoop.range]*fish.weights*params$dx)
																		
			
	##### COEFFICIENT OF VARIATION
	cv = sd(colMeans(fish.weights*slope.matrix[-zoop.range, c((0.5*params$N):(params$N-1))]*params					$dx))/
				mean(colMeans(fish.weights*slope.matrix[-zoop.range, c((0.5*params$N):(params$N-1))]*params					$dx))						
			
	##### ECOLOGICAL EFFICIENCY
	ecol.effic = (fish.biom/zoo.biom)/ee0
	
	sigma.store[i.mod,1] = fish.biom
	sigma.store[i.mod,2] = pb.store
	sigma.store[i.mod,3] = tp.store
	sigma.store[i.mod,4] = cv
	sigma.store[i.mod,5] = ecol.effic
	}
	
	if(i %in% ((l.m+l.s+1):N.runs)){
	i.mod <- (i-l.m-l.s)	
	params$sigma.z <- 1.75
	params$m <- 0
	params$K.z <- K.slide[(i-(l.m+l.s))]
	
	out <- PZF.solve(state = init.dist, parms = params, test = 1, zoo.beta = 1)
			
	slope.matrix <- matrix(unlist(out[1]), nrow = length(init.dist), ncol = params$N)
	feeding.zoo.matrix <- matrix(unlist(out[2]), nrow = length(params$x), ncol = params$N)
	feeding.fish.matrix <- matrix(unlist(out[3]), nrow = length(params$x), ncol = params$N)
	mort.zoo.matrix <- matrix(unlist(out[4]), nrow = length(params$x), ncol = params$N)
	mort.fish.matrix <- matrix(unlist(out[5]), nrow = length(params$x), ncol = params$N)		
		
	##### FISH BIOMASS
	slope.store = rowMeans(slope.matrix[,c((0.5*params$N):(params$N-1))])
	zoo.biom = sum(slope.store[zoop.range]*zoop.weights*params$dx)
	fish.biom = sum(slope.store[-zoop.range]*fish.weights*params$dx)
			
	##### P:B RATIO
	ave.mort.fish = rowMeans(mort.fish.matrix[(params$xfminref:params$xmaxref), 
																	c((0.5*params$N):(params$N-1))])
																	
	prod.store = sum(ave.mort.fish*slope.store[-zoop.range]*fish.weights*params$dx)

	pb.store = prod.store/fish.biom

	##### TOTAL THROUGHPUT
	ave.thru.fish = rowMeans(feeding.fish.matrix[(params$xfminref:params$xmaxref),																					 c((0.5*params$N):(params$N-1))])

	tp.store = sum(ave.thru.fish*slope.store[-zoop.range]*fish.weights*params$dx)
																		
			
	##### COEFFICIENT OF VARIATION
	cv = sd(colMeans(fish.weights*slope.matrix[-zoop.range, c((0.5*params$N):(params$N-1))]*params					$dx))/
				mean(colMeans(fish.weights*slope.matrix[-zoop.range, c((0.5*params$N):(params$N-1))]*params					$dx))						
			
	##### ECOLOGICAL EFFICIENCY
	ecol.effic = (fish.biom/zoo.biom)/ee0
	
	K.store[i.mod,1] = fish.biom
	K.store[i.mod,2] = pb.store
	K.store[i.mod,3] = tp.store
	K.store[i.mod,4] = cv
	K.store[i.mod,5] = ecol.effic	
	}
	

	
}

slope.store <- matrix(0, nrow = l.m*l.s*l.K, ncol = length(init.dist))
eval.store <- matrix(0, nrow = l.m*l.s*l.K, ncol = 1)
biom.store <- matrix(0, nrow = l.m*l.s*l.K, ncol = 3)
prod.store <- matrix(0, nrow = l.m*l.s*l.K, ncol = 3)
cv.store <- matrix(0, nrow = l.m*l.s*l.K, ncol = 3)
pb.store <- matrix(0, nrow = l.m*l.s*l.K, ncol = 3)
tp.store <- matrix(0, nrow = l.m*l.s*l.K, ncol = 3)


#pb <- txtProgressBar(min = 0, max = l.m, style = 3)

for(i in 1:l.m){
	params$m <- m.slide[i]
#	setTxtProgressBar(pb, i)
	
	for(j in 1:l.s){
		params$sigma.z <- sigma.slide[j]
		
		for(k in 1:l.K){
			curr.row = (i-1)*(l.s*l.K) + (j-1)*l.K + k # current row in storage matrices
			print(curr.row)
			if(curr.row < 1860 | curr.row > 3783){
			zoop.range = 1:(params$xzmaxref - params$xzminref + 1)
			zoop.weights = params$w[params$xzminref:params$xzmaxref]
			fish.weights = params$w[params$xfminref:params$xmaxref]
			params$K.z <- K.slide[k]
			out <- PZF.solve(state = init.dist, parms = params, test = 0, zoo.beta = 1)
			
			slope.matrix <- matrix(unlist(out[1]), nrow = length(init.dist), ncol = params$N)
			feeding.zoo.matrix <- matrix(unlist(out[2]), nrow = length(params$x), ncol = params$N)
			feeding.fish.matrix <- matrix(unlist(out[3]), nrow = length(params$x), ncol = params$N)
			mort.zoo.matrix <- matrix(unlist(out[4]), nrow = length(params$x), ncol = params$N)
			mort.fish.matrix <- matrix(unlist(out[5]), nrow = length(params$x), ncol = params$N)
			
			#### RESULTS STORAGE
			## Last year average size spectra slope
			slope.store[curr.row,] = rowMeans(slope.matrix[,c((0.5*params$N):(params$N-1))])
			
			## Largest e-val of slope.store
			#params$N = 2 # set number of time steps to 2 days for Newton-Raphson
			#eval.store[curr.row] = Newton(slope.store[curr.row,], z.b = 1)
			#params$N = 10*365 # set number of time steps back to 10 years	
			
			
			## Biomass density of fish and zooplankton
			biom.store[curr.row, 1] = sum(slope.store[curr.row, zoop.range]*zoop.weights*params$dx)
									# Average biomass density of zooplankton
			biom.store[curr.row, 2] = sum(slope.store[curr.row, -zoop.range]*fish.weights*params$dx)
									# Average biomass density of fish
			biom.store[curr.row, 3] = biom.store[curr.row, 1] + biom.store[curr.row, 2]
									
			## Coefficient of variation of fish and zooplankton
			cv.store[curr.row, 1] = 
				sd(colMeans(zoop.weights*slope.matrix[zoop.range, c((2*365):(params$N-1))]*params$dx))/
				mean(colMeans(zoop.weights*slope.matrix[zoop.range, c((2*365):(params$N-1))]*params$dx))										# Coefficient of variation of zooplankton biomass in final year
			cv.store[curr.row, 2] = 
				sd(colMeans(fish.weights*slope.matrix[-zoop.range, c((2*365):(params$N-1))]*params$dx))/
				mean(colMeans(fish.weights*slope.matrix[-zoop.range, c((2*365):(params$N-1))]*params$dx))										# Coefficient of variation of fish biomass in final year
			cv.store[curr.row, 3] = 
			sd(colMeans(c(zoop.weights,fish.weights)*slope.matrix[, c((2*365):(params$N-1))]*params$dx))/
		mean(colMeans(c(zoop.weights, fish.weights)*slope.matrix[, c((2*365):(params$N-1))]*params$dx))
					# Coefficient of variation of total system biomass in final year	
			
			## Production of fish and zooplankton
			ave.mort.zoo = rowMeans(mort.zoo.matrix[(params$xzminref:params$xzmaxref), 
																	c((2*365):(params$N-1))])
			ave.mort.fish = rowMeans(mort.fish.matrix[(params$xfminref:params$xmaxref), 
																	c((2*365):(params$N-1))])
			
			prod.store[curr.row, 1] = 																								sum(ave.mort.zoo*slope.store[curr.row,zoop.range]*zoop.weights*params$dx)
			prod.store[curr.row, 2] = 
							sum(ave.mort.fish*slope.store[curr.row, -zoop.range]*fish.weights*params$dx)
			prod.store[curr.row, 3] = prod.store[curr.row, 1] + prod.store[curr.row, 2]
			
			## Production to biomass ratio for fish and zooplankton
			pb.store[curr.row, 1] = prod.store[curr.row, 1]/biom.store[curr.row, 1]
			pb.store[curr.row, 2] = prod.store[curr.row, 2]/biom.store[curr.row, 2]
			pb.store[curr.row, 3] = pb.store[curr.row, 1] + pb.store[curr.row, 2]
			
			# Throughput of fish and zooplankton
			ave.thru.zoo = rowMeans(feeding.zoo.matrix[(params$xzminref:params$xzmaxref),																					 c((2*365):(params$N-1))])
			ave.thru.fish = rowMeans(feeding.fish.matrix[(params$xfminref:params$xmaxref),																					 c((2*365):(params$N-1))])
			tp.store[curr.row, 1] = sum(ave.thru.zoo*slope.store[curr.row, zoop.range]*
																		zoop.weights*params$dx)
			tp.store[curr.row, 2] = sum(ave.thru.fish*slope.store[curr.row, -zoop.range]*
																		fish.weights*params$dx)
			tp.store[curr.row, 3] = tp.store[curr.row, 1] + tp.store[curr.row, 2]
			}
		
	}
}

}
 
write.matrix(slope.store, file = "slope.csv", sep = ",")
write.matrix(biom.store, file = "biom.csv", sep = ",")
write.matrix(eval.store, file = "eval.csv", sep = ",")
write.matrix(cv.store, file = "cv.csv", sep = ",")
write.matrix(prod.store, file = "prod.csv", sep = ",")
write.matrix(pb.store, file = "pb.csv", sep = ",")
write.matrix(tp.store, file = "tp.csv", sep = ",") 