library(raster)
library(sp)
library(rgeos)

# Read in Hardwood plots # 
HWplots <- dget("Data/Hardwood_plots.txt")

###########################################################
# Save only count data that falls within Hardwood plots ### 
# spec.mat <- spec.mat[HWplots,,,]
# date <- date[HWplots,,]
# time <- time[HWplots,,]
#
# Save these new data sets
# saveRDS(spec.mat,"Data/Abundance_array.rds")
# saveRDS(time,"Data/TimeOfCount.rds")
# saveRDS(date,"Data/DateOfCount.rds")
#
###########################################################

# Read in pre-processed data 

# Abundance data
spec.mat <- readRDS("Data/Abundance_array.rds")

# Time of count (standardized) 
time <- readRDS("Data/TimeOfCount.rds")

# Date of count (standardized) 
date <- readRDS("Data/DateOfCount.rds")

# Mean Temperatures 
mayMeanTemp <- readRDS("Data/mayMeanTemp.rds")
juneMeanTemp <- readRDS("Data/juneMeanTemp.rds")

juneMeanTemp <- (juneMeanTemp - mean(juneMeanTemp,na.rm = TRUE))/sd(juneMeanTemp,na.rm = TRUE)

# Precip
precip <- readRDS("Data/precip.rds")

# collapse months into years by sum precip #
Precip <- apply(precip,c(1,2),sum)

# Standardize 
Precip <- (Precip-mean(Precip,na.rm = TRUE))/sd(Precip,na.rm = TRUE)

# Thermal Sum 
ThermSum <- readRDS("Data/thermalsum.rds")

# Structural Complexity 
complex <- raster::shapefile("Spatial_Layers/SurveyPoints/Schwarz_CanopyComplex.shp")
CanopyComplex <- complex$CnpyCmp
# standardize canopy 
StrComplex <- (CanopyComplex - mean(CanopyComplex,na.rm = TRUE))/sd(CanopyComplex, na.rm = TRUE)
StrComplex[is.na(StrComplex)]<-0

############################################################################################################
#
# Spatial Proximity 
# 
############################################################################################################

library(raster)
library(spdep)

HWpoints <- raster::shapefile("Spatial_Layers/SurveyPoints/HWsurveypoints.shp")

# Compute the neighborhood (2nd order)
neigh <- dnearneigh(HWpoints,d1 = 0, d2 = 501)

plot(neigh,HWpoints@coords)

# Store the neighborhood in the format needed by the BUGS model
cardneigh <- card(neigh)

adjmat <- matrix(nrow = length(neigh), ncol = max(cardneigh))
for (i in 1:length(neigh)) {
    adjmat[i, 1:cardneigh[i]] <- neigh[[i]]
}

###########################################################################################################
#
# Save data for analysis 
#
###########################################################################################################
# generate a object that defines the unique combinations of species # 
upperd <- which(upper.tri(array(NA,c(13,13)))==TRUE,arr.ind = TRUE)

# Abundance # 
win.data <- list(y = spec.mat,
                 nschwarz = dim(spec.mat)[1],
                 nyears = dim(spec.mat)[3],
                 nreps = dim(spec.mat)[2],
                 nspp = dim(spec.mat)[4],
                 spp = 1:13,
                 ncovs = 4,
                 pcovs = 2,
                 ThermalSum = ThermSum,
                 StrComplex = StrComplex,
                 Precip = Precip, 
                 JuneTemp = juneMeanTemp,
                 cardneigh = cardneigh,
                 adjmat = adjmat,
                 date = date,
                 time = time,
                 upperd = upperd)

# saveRDS(win.data,"Data/AbunData.rds")

# CONVERT ABUNDANCE DATA TO OCCUPANCY #
occ.data <- win.data
occ.data$y[occ.data$y>1]<-1        

# saveRDS(occ.data,"Data/OccuData.rds")      