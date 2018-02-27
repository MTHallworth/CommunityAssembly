library(raster)
library(jagsUI)

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

params<-c("gam","ranked.order","ranked.vw","N",
          "beta","betaP","MuBeta","MuGam","meanDiffBeta","meanDiffGams")

n.iter <- 10
n.burn <- 5
n.thin <- 1
n.adapt <- 0
n.chains <- 3

 Sys.time()
 a<-proc.time()
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
proc.time()-a
 Sys.time()

saveRDS(dm.fit.spp,"dm.fit.spp.rds")

source("Scripts/sims.list.R")

Sys.time()
a<-proc.time()
dm <- process.output(dm.fit.spp[[1]],Rhat = FALSE)
proc.time()-a
Sys.time()