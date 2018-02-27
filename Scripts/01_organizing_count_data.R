
# Read in valley-wide data set

VWdata <- read.csv("Data/Valley-wide_bird_census_1999_2016.csv", header = TRUE, sep = ",")


#Let's look at the first few rows of data 

head(VWdata)

# Remove data that needs to be deleted # This was fixed in FileMaker #

VWdata<-VWdata[VWdata$Comments!="**PLEASE DELETE**",]

# Some of the species codes have new lines for example "BTBW\n" this next code removes that & turns "BTBW\n" into "BTBW"

VWdata$Species<-gsub(pattern = "[\r\n]", replacement = "", x = VWdata$Species)
VWdata$Species<-gsub(pattern = "[\r\v]", replacement = "", x = VWdata$Species)

# Some rows of the data were entered by accident with no new species - they have blank species#
# remove those rows from the data set #
VWdata<-VWdata[VWdata$Species != "", ]

# remove data that were used for a separate analysis - Multiple Observer experiment #
VWdata <- VWdata[VWdata$Replicate != "4-MultObs",]
VWdata$Replicate[VWdata$Replicate == "-3"] <- 3
VWdata$Replicate <- droplevels(VWdata$Replicate)

# change replicate to numeric
VWdata$Replicate <- as.numeric(as.character(VWdata$Replicate))

# Use only the data collected within 50m of the point so no double counting #
VWdata<-VWdata[VWdata$Distance==1,]

# Determine the Number of counts
replicate <-min(VWdata$Replicate,na.rm=TRUE):max(VWdata$Replicate,na.rm=TRUE)

# Determine the number of periods within a count # - should be 3 - one record has value 4
# Identify and fix the record that was entered incorrectly - NOTE this was updated in FileMaker by MTH #
 
period<-min(VWdata$Period,na.rm=TRUE):max(VWdata$Period,na.rm=TRUE) # Within the count

# number of schwarz plots surveyed #
sites<-levels(as.factor(VWdata$Plot))

# number of species observed #
spp <- levels(factor(VWdata$Species))

# List the species of interest (S.O.I) #
S.O.I <- c("MYWA","BTHW","BHVI","BLPW","BTBW","MAWA","OVEN","BLBW","AMRE","REVI","CAWA","BAWW","NAWA")

# put them in ABC order
S.O.I <- S.O.I[order(S.O.I)]

# Format the dates
VWdata$Date2 <- as.Date(VWdata$Date, format="%m/%d/%Y")

# observation year #
VWdata$Year <- as.numeric(format(VWdata$Date2,"%Y"))

# vector of years
years <- 1999:2015

# observation ordinal day #
VWdata$DOY <- as.numeric(format(VWdata$Date2, "%j"))

# Create a Plot x Replicate array of abundance for the species of interest #
# this makes an empty 3 dimensional array - Plot x Replicate x Species that is filled with NAs #

# create an empty array (filled with NA) -   # of rows = sites   # of columns = periods   # of dimensions = replicates
spec.mat <- array(NA, c(length(sites),4,length(years),length(S.O.I)))

rownames(spec.mat)<-sites
colnames(spec.mat) <- replicate[1:4]
dimnames(spec.mat)[3][[1]] <- years
dimnames(spec.mat)[4][[1]] <- S.O.I

# Here is how to fill the array with the data that we want # 
for(y in 1:length(years)){
temp.yr <- subset(VWdata, Year == years[y])
for(s in 1:length(S.O.I)){
temp.sp <- subset(temp.yr, Species == S.O.I[s]) 
# for each replicate #
  for(i in 1:4){
      temp.rp <- subset(temp.sp, Replicate == replicate[i])
# for each plot #
         for(k in sites){
               temp.st <- subset(temp.rp, Plot == as.numeric(k))
               if(years[y]>2005){
               spec.mat[k,i,y,s] <- as.numeric(length(temp.st$Species[temp.st$New.Record == 1]))}else{spec.mat[k,i,y,s] <- as.numeric(length(temp.st$Species))}
          }  # For plot
      }  # For Replicate
  }# For Species
}# years 

### OBSERVATION COVARIATES 

# DATE OF COUNT

doy <- array(NA, c(length(sites), length(replicate), length(years)))

rownames(doy) <- sites
dimnames(doy)[2][[1]] <- replicate
dimnames(doy)[3][[1]] <- years

for(y in 1:length(years)){
temp.yr <- subset(VWdata, Year == years[y])
for(k in sites){
  for(i in 1:4) {
    temp <- subset(temp.yr, Replicate == replicate[i] & Plot==k)
    temp.DOY <- as.numeric(levels(factor(temp$DOY)))
    if (length(temp.DOY)>0){
        doy[k,i,y] <- min(temp.DOY)
        } # if
      } # i 
   } # k
} # y

# Change the array from characters to numeric #
doy <- structure(as.numeric(doy), dim=dim(doy), dimnames=dimnames(doy))

# All covariates need to be standardized. You can use the function scale() to do that #
doy.scaled<-scale(doy)
doy.scaled[is.na(doy.scaled)]<-0

# FORMAT TIME OF COUNTS #

VWdata$Time2 <-gsub(" AM", "", VWdata$Time)

# Create a site x occasion matrix of times 

time <- array(NA, c(length(sites), length(replicate), length(years)))

rownames(time) <- sites
colnames(time)<-replicate
dimnames(time)[3][[1]] <- years

for(y in 1:length(years)){
temp.yr <- subset(VWdata,Year == years[y])
for(k in sites){
for(i in 1:4){
    temp <- subset(temp.yr, Replicate== replicate[i] & Plot==k)
    temp.time <- temp$Time2
    if(length(temp.time) > 0) {
        temp.time <- substr(temp.time,1,5)
        time[k,i,y] <- as.numeric(gsub(":","",(unique(temp.time)[1])))
        }
    }
}
}

# Scale the time covariate #
time <- (time - mean(time,na.rm = TRUE))/sd(time,na.rm = TRUE)
# Fill missing data with the mean - which is 0 after being scaled #
time[is.na(time)]<-0


# 2011 was a little weird counts were numbered 0,1,2 - instead of 1,2,3 
# We need to shift those values over to match other years #
spec.mat[,2:4,13,]<-spec.mat[,1:3,13,]

# Now save only the counts that were labeled 1,2,3 in all years #
spec.mat <- spec.mat[,2:4,,]
date <- doy[,2:4,]
time <- time[,2:4,]

# standardize date and time # 
date <- (date - mean(date,na.rm = TRUE))/sd(date,na.rm = TRUE)
date[is.na(date)] <- 0

time <- (time - mean(time,na.rm = TRUE))/sd(time,na.rm = TRUE)
time[is.na(time)] <- 0



# SAVE THE DATA TO FILE # 
str(spec.mat, 1)
str(time,1)
str(date,1)
#saveRDS(spec.mat,"Data/Abundance_array.rds")
#saveRDS(time,"Data/TimeOfCount.rds")
#saveRDS(date,"Data/DateOfCount.rds")

