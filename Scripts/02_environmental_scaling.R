library(raster)
library(sp)
library(prism)
library(maptools)


boundary<-shapefile("Spatial_Layers/HBEFboundary/HBEFboundary_WGS84.shp")
elev<-raster("Spatial_Layers/DEM/hb10mdem.txt")
slope<-terrain(elev,"slope")
aspect<-terrain(elev,"aspect")
roughness <- terrain(elev,"roughness")
HWpoints <- shapefile("Spatial_Layers/SurveyPoints/HWsurveypoints.shp")

#############################################################################################################
# 
# Statistical downscaling of PRISM data from 800m to 50m 
#
#############################################################################################################
#
#  TEMPERATURE DATA 
#
#############################################################################################################
#
#
# HERE YOU NEED PRISM DATA THAT ARE ON EXTERNAL HARD DRIVE #
# 
daily <- raster("E:/PRISMdata/daily2014/tmax/PRISM_tmax_stable_4kmD1_20140601_bil.bil")
monthly<-raster("E:/PRISMdata/PRISM_tmax_30yr_normal_800mM2_06_bil/PRISM_tmax_30yr_normal_800mM2_06_bil.bil")

# project data 

boundary<-spTransform(boundary,CRS=crs(daily))
daily<-crop(daily,boundary)
monthly<-crop(monthly,boundary)

UTM19N<-"+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

daily<-projectRaster(daily,crs=UTM19N)
monthly<-projectRaster(monthly,crs=UTM19N)

boundary<-spTransform(boundary,CRS=UTM19N)
elev<-projectRaster(elev,crs=UTM19N)

# Crop to the HBEF boundary

daily<-crop(daily,boundary)
monthly<-crop(monthly,boundary)

# Extract elevation 
elevDaily<-rasterToPolygons(daily)
elevDaily$Elev<-extract(elev,elevDaily,mean,na.rm=TRUE)

# read in the 50x50m grid of HBEF
HBEFgridcovs<-read.csv("Spatial_Layers/HBEFgridcovs.csv")
HBEFgridcovs_sp<-SpatialPoints(cbind(HBEFgridcovs[,1:2]))

# Set crs
crs(HBEFgridcovs_sp)<-UTM19N

# Extract data to grid #
monthly.grid<-extract(monthly,HBEFgridcovs_sp)
daily.grid<-extract(daily,HBEFgridcovs_sp)
elev.grid<-extract(elev,HBEFgridcovs_sp,buffer=50,fun=mean)
Zelev.grid<-extract(elevDaily,HBEFgridcovs_sp)$Elev

############################################################################################################
# The number of grid cells 
n.grids<-20844

# Vector of years
yrs <- 1998:2015

# Generate empty list
May <- vector('list',length(yrs))

# Start loop for years # 
for(y in 1:length(yrs)){

# read in files that match year #
min_bils<-list.files(pattern="*bil.bil",paste0("E:/PRISMdata/MeanMayTemp_",yrs[y]),full.names=T)
min_bils<-min_bils[!min_bils %in% list.files(pattern=".aux.xml",paste0("E:/PRISMdata/MeanMayTemp_",yrs[y]),full.names=T)]

# stack the rasters #
layers_May<-stack(sapply(min_bils,raster))

# Crop to boundary #
b<-spTransform(boundary,CRS=crs(layers_May))
stack_hbef_may<-crop(layers_May,b)
stack_hbef_may<-projectRaster(stack_hbef_may,crs=UTM19N)

# for each day get the values then downscale to grid point #
days1 <- nlayers(layers_May)
tempsMay <- dailyMay <- array(NA,c(nrow(HBEFgridcovs),days1))

for(i in 1:days1){
  tempsMay[,i]<-extract(stack_hbef_may[[i]],HBEFgridcovs_sp)
  dailyMay[,i]<-tempsMay[,i]+(Zelev.grid-elev.grid)*0.006554
}

# save into year list 
May[[y]] <- dailyMay
}

MayMean <- vector("list",length(yrs))
for(y in 1:length(yrs)){
MayMean[[y]] <- raster::rasterFromXYZ(cbind(HBEFgridcovs_sp@coords,apply(May[[y]],1,mean,na.rm = TRUE)))
}

June <- vector('list',length(yrs))
for(y in 1:length(yrs)){
min_bils<-list.files(pattern="*bil.bil",paste0("E:/PRISMdata/MeanJuneTemp_",yrs[y]),full.names=T)
min_bils<-min_bils[!min_bils %in% list.files(pattern=".aux.xml",paste0("E:/PRISMdata/MeanJuneTemp_",yrs[y]),full.names=T)]

layers_June<-stack(sapply(min_bils,raster))
b<-spTransform(boundary,CRS=crs(layers_June))
stack_hbef_june<-crop(layers_June,b)
stack_hbef_june<-projectRaster(stack_hbef_june,crs=UTM19N)

days1 <- nlayers(layers_June)
tempsJune <- dailyJune <- array(NA,c(nrow(HBEFgridcovs),days1))

for(i in 1:days1){
  tempsJune[,i]<-extract(stack_hbef_june[[i]],HBEFgridcovs_sp)
  dailyJune[,i]<-tempsJune[,i]+(Zelev.grid-elev.grid)*0.006554
}

June[[y]] <- dailyJune
}

JuneMean <- vector("list",length(yrs))
for(y in 1:length(yrs)){
JuneMean[[y]] <- raster::rasterFromXYZ(cbind(HBEFgridcovs_sp@coords,apply(June[[y]],1,mean,na.rm = TRUE)))
}

#saveRDS(stack(JuneMean),"Data/juneMeanTemp.rds")
#saveRDS(stack(MayMean),"Data/mayMeanTemp.rds")

#############################################################################################################
# 
# Statistical downscaling of PRISM data from 800m to 50m 
#
#############################################################################################################
#
#  Precipitation DATA 
#
#############################################################################################################

#scaling factor in HBEF of elevation and precipitation # 
# read precipitation data #

psd <- read.csv("Data/psd_.txt")

# change no.data value to NA #
psd[psd == -99] <- NA

psd$date <- as.POSIXlt(psd[,1],"%Y-%M-%D")
psd$month <- format(psd$date,"%b")
psd$year <- as.numeric(format(psd$date,"%Y"))

# Read in rain gage locations to extract elevation
rg <- shapefile("Spatial_Layers/HBEFraingage/hbef_raingage.shp")

rg <- spTransform(rg,CRS(UTM19N))

rg_elev <- data.frame(rgID = rg$ID,
                      rgElev = raster::extract(elev,rg@coords[,1:2]))

rg_elev$IDnum <- as.numeric(gsub("RG","",rg_elev$rgID))

rg_elev <- rg_elev[order(rg_elev$IDnum),]

rainday <- psd[apply(psd[,c(2:18,20:26)],1,sum,na.rm = TRUE)>0,]

rainday <- subset(rainday, subset =(month %in% c("May","Jun","Jul") & year > 1998))

pptBeta <- apply(rainday[,c(2:18,20:26)],1,FUN = function(x){coef(lm(x~rg_elev$rgElev))[2]})

pptDaily <- raster ("E:/DailyPrecip/PRISM_ppt_30yr_normal_4kmM2_05_asc/PRISM_ppt_30yr_normal_4kmM2_05_asc.asc")
pptMonthly<-raster("E:/DailyPrecip/PRISM_ppt_30yr_normal_800mM2_05_asc/PRISM_ppt_30yr_normal_800mM2_05_asc.asc")

boundary<-spTransform(boundary,CRS=crs(pptDaily))
pptDaily<-crop(pptDaily,boundary)
pptMonthly<-crop(pptMonthly,boundary)

UTM19N<-"+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

pptDaily<-projectRaster(pptDaily,crs=UTM19N)
pptMonthly<-projectRaster(pptMonthly,crs=UTM19N)

boundary<-spTransform(boundary,CRS=UTM19N)
elev<-projectRaster(elev,crs=UTM19N)

elevDaily<-rasterToPolygons(pptDaily)
#elevDaily <- rgeos::gIntersection(elevDaily,boundary,byid = TRUE)
elevDaily$Elev<-extract(elev,elevDaily,mean,na.rm=TRUE)

HBEFgridcovs<-read.csv("Spatial_Layers/HBEFgridcovs.csv")
HBEFgridcovs_sp<-SpatialPoints(cbind(HBEFgridcovs[,1:2]))

crs(HBEFgridcovs_sp) <- crs(boundary)

monthly.grid<-extract(pptMonthly,HBEFgridcovs_sp)
daily.grid<-extract(pptDaily,HBEFgridcovs_sp)
elev.grid<-extract(elev,HBEFgridcovs_sp,buffer=50,fun=mean)
Zelev.grid<-extract(elevDaily,HBEFgridcovs_sp)$Elev

years <- 1998:2015

mons <- c("05","06","07")
precipMon <- scaledPrecip <- array(NA,c(20844,length(years),length(mons)))
precip <- array(NA,c(242,length(years),length(mons)))

HWplots <- rgeos::gBuffer(HWpoints, width = 50, byid = TRUE)

for(y in 1:length(years)){
for(m in 1:length(mons)){
precipRast <- raster(paste0("E:/Precipitation/PRISM_ppt_stable_4kmM3_",years[y],mons[m],"_bil.bil"))
precipRast <- projectRaster(precipRast,crs = UTM19N)
  precipMon[,y,m]<-extract(precipRast,HBEFgridcovs_sp)
  scaledPrecip[,y,m]<-precipMon[,y,m]+(Zelev.grid-elev.grid)*mean(pptBeta)

tempRast <- rasterFromXYZ(cbind(HBEFgridcovs_sp@coords,scaledPrecip[,y,m]))
precip[,y,m] <- extract(tempRast,HWplots,mean)
}
}

#saveRDS(precip,"Data/precip.rds")
########################################################################################

# Thermal Sum # 
ThermalSum <- readRDS("C:/Users/hallworthm/Dropbox (Smithsonian)/Cost_of_Reproduction/Spatial_Layers/ThermalSumYr.rds")
str(ThermalSum)
June1_TS <- apply(ThermalSum[,1:153,4:21],c(1,3),sum)

ThermalSum <- array(NA,c(242,length(years)))
for(y in 1:length(years)){
tempRast <- rasterFromXYZ(cbind(HBEFgridcovs_sp@coords,June1_TS[,y]))
ThermalSum[,y] <- extract(tempRast,HWplots,mean)
}

#saveRDS(ThermalSum,"Data/thermalsum.rds")
