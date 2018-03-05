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


# Read in abundance data #
win.data <- readRDS("Data/AbunData.rds")
win.data$pi <- pi

# Function to set up the initial values
## ------------------------------------------------------------------------
N<-function(x){
  Nst <- array(NA,dim(x)[c(1,3:4)])
  nspp <- dim(x)[4]
  for(s in 1:nspp){
    Nst[,,s]<-apply(x[,,,s],c(1,3),max,na.rm=TRUE)+3
  }
  
  # If Nst==NA change value to 3 #
  Nst[Nst==-Inf]<-NA
  Nst[is.na(Nst)]<-3
  return(Nst)
}

## ---- tidy = FALSE, eval = FALSE-----------------------------------------
 
inits<-function() list(N=N(win.data$y))

params<-c("MuAlpha","MupInt","MuError",
          "MuBeta","sigmaBeta","beta","MuBetaP",
          "sigmaBetaP","betaP","MuGam","sigmaGam",
          "gam", "pInt","alpha","Error","phi.s",
          "ranked.order","ranked.vw","N","recruits.ha",
          "Ntot","log.abund","meanDiffBeta","meanDiffGams")

n.iter <- 50000
n.burn <- 10000
n.thin <- 10
n.adapt <- 1000
n.chains <- 3

 Sys.time()
 a<-Sys.time()
 dm.fit.spp<-jagsUI::jags.basic(model ="Models/dm_eco_spatial_alt.txt", 
                                module = c('glm'),
								data = win.data,
								inits = inits,
								parameters.to.save = params,
								parallel = TRUE,
								n.adapt = n.adapt,
								n.iter = n.iter,
								n.burnin = n.burn,
								n.thin = n.thin,
								n.chains = n.chains,
								save.model = TRUE,
								DIC = FALSE,
			                    seed = 04823)
Sys.time()-a
Sys.time()

saveRDS(dm.fit.spp,"dm_fit_spp_basic.rds")
				
source("Scripts/sims.list.R")	

Sys.time()
a<-Sys.time()
dm <- process.output(dm.fit.spp[[1]],Rhat = FALSE)
Sys.time()-a
Sys.time()

# save results #
saveRDS(dm[c(2:12)],paste0("DynAbun_eco_spatial_",Sys.Date(),".rds"))
saveRDS(dm[1],paste0("DynAbun_simslist_",Sys.Date(),".rds"))

