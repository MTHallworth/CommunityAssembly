library(raster)
library(jagsUI)

# Read in abundance data #
win.data <- readRDS("Data/OccuData.rds")
win.data$pi <- pi


# Inits function for occupancy model #
Zst <- array(NA,c(win.data$nschwarz,win.data$nyears,win.data$nspp))
for(s in 1:win.data$nspp){
Zst[,,s] <- apply(win.data$y[,,,s],c(1,3),max,na.rm = TRUE)
}
Zst[Zst == "-Inf"]<-1

inits <- function(){list(z = Zst)}


params <- c("m.dispersal","v.dispersal","alpha.lgamma", "beta.lgamma",
             "n.occupied","growthrate","turnover","spp.rich", "mean.dispersal",
             "Yrs.Occ","beta.psi","alpha.psi")


Sys.time()
a<-proc.time()
Dyn.Occ <- jagsUI::jags.basic(data = win.data,
							model = "Models/dynocc_eco_spatial_alt.txt",
							parameters.to.save = params,
							inits = inits,
							n.chains = 3,
							n.iter = 10,
							n.adapt = 0,
							n.thin = 1, 
							n.burnin = 5,
							save.model = TRUE,
							parallel = TRUE,
							seed = 04823)
													
					
source("Scripts/sims.list.R")	
	

DynOcc <- process.output(Dyn.Occ[[1]],Rhat = FALSE)
proc.time()-a
Sys.time()