#### -------------------- set path to library dir --------------------------- #### 

library(parallel,lib.loc="/home/hallworthm/CostOfRepro/Packages")
library(coda, lib.loc="/home/hallworthm/CostOfRepro/Packages")
library(grid)
library(lattice,lib.loc="/home/hallworthm/CostOfRepro/Packages")
library(rjags,lib.loc="/home/hallworthm/CostOfRepro/Packages")
library(jagsUI,lib.loc="/home/hallworthm/CostOfRepro/Packages")


#### -------------------- set path to library dir --------------------------- ####

my.path <- "/home/hallworthm/CostOfRepro/Packages"

.libPaths(c(my.path,.libPaths()))


# Read in data #

# Read in abundance data #
win.data <- readRDS("Data/OccuData.rds")
win.data$pi <- pi

# Inits function for occupancy model #
Zst <- array(NA,c(win.data$nschwarz,win.data$nyears,win.data$nspp))
for(s in 1:win.data$nspp){
Zst[,,s] <- apply(win.data$occ[,,,s],c(1,3),max,na.rm = TRUE)
}
Zst[Zst == "-Inf"]<-1

inits <- function(){list(z = Zst)}


params <- c(# Hyperparameters 
            "Mulgamma","sigmalgamma",
            "MuBetaPsi","sigmaBetaPsi",
            "MuBetaP", "sigmaBetaP",
            "MuAlphaPsi", "MuAlphaP",
            "MuBetalphi", "MuAlphaLgamma",
            "sigmaAlphaPsi","sigmaAlphaP",
            "sigmaBetalphi", "sigmaAlphaLgamma",
            "sigmaAlphaLphi", 
            # estimates 
            "alpha.psi", "alpha.p","beta.lphi","beta.lgamma",
            "beta.psi", "beta.p","phi","alpha.lgamma",
            "m.dispersal","v.dispersal","mean.dispersal",
             "n.occupied","Yrs.Occ","growthrate","turnover","spp.rich", 
            "effectlgamma","meanDiffBeta")




n.iter <- 50000
n.burn <- 10000
n.thin <- 10
n.adapt <- 5000
n.chains <- 3

Sys.time()
a<-Sys.time()
Dyn.Occ <- jagsUI::jags.basic(data = win.data,
							model = "Models/dynocc_eco_spatial_alt.txt",
							parameters.to.save = params,
							inits = inits,
							n.chains = n.chains,
							n.iter = n.iter,
							n.adapt = n.adapt,
							n.thin = n.thin, 
							n.burnin = n.burn,
							save.model = TRUE,
							parallel = TRUE,
							seed = 04823)
													
Sys.time()-a
Sys.time()	

# Save model output incase need to re-run or update #
saveRDS(Dyn.Occ,"dynocc_basic.rds")

a<-Sys.time()				
source("Scripts/sims.list.R")	
	
DynOcc <- process.output(Dyn.Occ[[1]],Rhat = FALSE)

Sys.time()-a
Sys.time()

# save results #
saveRDS(DynOcc[c(2:12)],paste0("DynOcc_eco_spatial_",Sys.Date(),".rds"))
saveRDS(DynOcc[1],paste0("DynOcc_simslist_",Sys.Date(),".rds"))


