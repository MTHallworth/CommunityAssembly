#########################################################################
#
#
#                ABUNDANCE POPULATION LEVEL METRIC 
#
#
#########################################################################

# The input for abundance needs to be called dm - needs to be sims.list # 

betaOrder <- order(c(dm$mean$MuBeta,dm$mean$MuGam[2:3]))

# Combine relevant data into a single array #
betaSL <- cbind(dm$sims.list$MuBeta,dm$sims.list$MuGam[,2:3])
colnames(betaSL) <- c("ThermSum","StrComplex","Precip","JuneTemp","Density","ConSpecific")

paramOrd <- rev(order(apply(betaSL,2,median)))

#tiff("Figures/Abundance_example.tiff",res = 600, width = 3600, height = 3600, units = "px", compress = "lzw")
# Joy Division Plot #
par(mar = c(4,6,3,3), bty = "l")
plot(NA,type = "n", 
     ylim = c(1,7.5),
     xlim = c(-1,1),
     xlab = "",
     ylab = "",
     yaxt = "n",
     xaxt = "n")
axis(1, cex.axis = 0.75)
axis(2, las = 2, at = 1:6, labels = colnames(betaSL)[paramOrd], cex = 0.8)
mtext(side = 1, at = NA, line = 2, cex = 0.75,
      text = expression(paste("Posterior ",beta," estimate")))

scalepar <- 0.2

for(i in 6:1){

x <- density(betaSL[,paramOrd[i]])$x
y <- density(betaSL[,paramOrd[i]])$y*scalepar

# FULL DENSITY POLYGON #
polygon(x = x,
        y = y+i,
        col = "gray80",
        border = "gray88")


# 95% DENSITY POLYGON #
polygon( x = c(x[findInterval(quantile(betaSL[,paramOrd[i]],probs = 0.05),x)],
               x[c(findInterval(quantile(betaSL[,paramOrd[i]],probs = 0.05),x):
                   findInterval(quantile(betaSL[,paramOrd[i]],probs = 0.95),x))],
               x[findInterval(quantile(betaSL[,paramOrd[i]],probs = 0.95),x)]),
         y = c(i,(y[c(findInterval(quantile(betaSL[,paramOrd[i]],probs = 0.05),x):
                 findInterval(quantile(betaSL[,paramOrd[i]],probs = 0.95),x))]+i),i),
         col = "gray60",
         border = "gray88")
# 
segments(x0 = median(betaSL[,paramOrd[i]]),
         y0 = max(y)+i,
         x1 = median(betaSL[,paramOrd[i]]),
         y1 = i, 
         col = "white")
}
abline(v = 0, col = "gray88",lwd = 1.5,lty = 3)
#dev.off()
system("open Figures//Abundance_example.tiff")




#########################################################################
#
#
#                OCCUPANCY POPULATION LEVEL METRIC 
#
#
#########################################################################

DYN <- readRDS("E:/DynOcc_simslist_2018-03-06.rds")

betaOrder <- order(apply(DYN$sims.list$Mulgamma,2,median))

# Combine relevant data into a single array #
paramnames <- c("Density","ThermSum","StrComplex","Precip","JuneTemp")

#tiff("Figures/Occupancy_populationMetrics.tiff",res = 600, width = 3600, height = 3600, units = "px", compress = "lzw")
# Joy Division Plot #
par(mar = c(4,6,3,3), bty = "l")
plot(NA,type = "n", 
     ylim = c(1,6.5),
     xlim = c(-0.5,5),
     xlab = "",
     ylab = "",
     yaxt = "n",
     xaxt = "n")
axis(1, cex.axis = 0.75)
axis(2, las = 2, at = 1:5, labels = paramnames[betaOrder], cex = 0.8)
mtext(side = 1, at = NA, line = 2, cex = 0.75,
      text = expression(paste("Posterior ",beta," estimate")))

scalepar <- 0.2

for(i in 5:1){

x <- density(DYN$sims.list$Mulgamma[,betaOrder[i]])$x
y <- density(DYN$sims.list$Mulgamma[,betaOrder[i]])$y*scalepar

# FULL DENSITY POLYGON #
polygon(x = x,
        y = y+i,
        col = "gray80",
        border = "gray88")


# 95% DENSITY POLYGON #
polygon( x = c(x[findInterval(quantile(DYN$sims.list$Mulgamma[,betaOrder[i]],probs = 0.05),x)],
               x[c(findInterval(quantile(DYN$sims.list$Mulgamma[,betaOrder[i]],probs = 0.05),x):
                   findInterval(quantile(DYN$sims.list$Mulgamma[,betaOrder[i]],probs = 0.95),x))],
               x[findInterval(quantile(DYN$sims.list$Mulgamma[,betaOrder[i]],probs = 0.95),x)]),
         y = c(i,(y[c(findInterval(quantile(DYN$sims.list$Mulgamma[,betaOrder[i]],probs = 0.05),x):
                 findInterval(quantile(DYN$sims.list$Mulgamma[,betaOrder[i]],probs = 0.95),x))]+i),i),
         col = "gray60",
         border = "gray88")
# 
segments(x0 = median(DYN$sims.list$Mulgamma[,betaOrder[i]]),
         y0 = max(y)+i,
         x1 = median(DYN$sims.list$Mulgamma[,betaOrder[i]]),
         y1 = i, 
         col = "white")
}
abline(v = 0, col = "gray88",lwd = 1.5,lty = 3)
#dev.off()
system("open Figures/Occupancy_populationMetrics.tiff")



###########################################################################################################
#
# DIFFERENCE BETWEEN SPECIES 
#
###########################################################################################################
betaOrder <- order(apply(DYN$sims.list$meanDiffBeta,2,median))

# Combine relevant data into a single array #
paramnames <- c("Density","ThermSum","StrComplex","Precip","JuneTemp")

tiff("Figures/Occupancy_populationMetrics_niche.tiff",res = 600, width = 3600, height = 3600, units = "px", compress = "lzw")
# Joy Division Plot #
par(mar = c(4,6,3,3), bty = "l")
plot(NA,type = "n", 
     ylim = c(1,6.5),
     xlim = c(-0.5,1.5),
     xlab = "",
     ylab = "",
     yaxt = "n",
     xaxt = "n")
axis(1, cex.axis = 0.75)
axis(2, las = 2, at = 1:5, labels = paramnames[betaOrder], cex = 0.8)
mtext(side = 1, at = NA, line = 2, cex = 0.75,
      text = expression(paste("Posterior ",beta," estimate")))

scalepar <- 0.2

for(i in 5:1){

x <- density(DYN$sims.list$meanDiffBeta[,betaOrder[i]])$x
y <- density(DYN$sims.list$meanDiffBeta[,betaOrder[i]])$y*scalepar

# FULL DENSITY POLYGON #
polygon(x = x,
        y = y+i,
        col = "gray80",
        border = "gray88")


# 95% DENSITY POLYGON #
polygon( x = c(x[findInterval(quantile(DYN$sims.list$meanDiffBeta[,betaOrder[i]],probs = 0.05),x)],
               x[c(findInterval(quantile(DYN$sims.list$meanDiffBeta[,betaOrder[i]],probs = 0.05),x):
                   findInterval(quantile(DYN$sims.list$meanDiffBeta[,betaOrder[i]],probs = 0.95),x))],
               x[findInterval(quantile(DYN$sims.list$meanDiffBeta[,betaOrder[i]],probs = 0.95),x)]),
         y = c(i,(y[c(findInterval(quantile(DYN$sims.list$meanDiffBeta[,betaOrder[i]],probs = 0.05),x):
                 findInterval(quantile(DYN$sims.list$meanDiffBeta[,betaOrder[i]],probs = 0.95),x))]+i),i),
         col = "gray60",
         border = "gray88")
# 
segments(x0 = median(DYN$sims.list$meanDiffBeta[,betaOrder[i]]),
         y0 = max(y)+i,
         x1 = median(DYN$sims.list$meanDiffBeta[,betaOrder[i]]),
         y1 = i, 
         col = "white")
}
abline(v = 0, col = "gray88",lwd = 1.5,lty = 3)
#dev.off()
system("open Figures/Occupancy_populationMetrics_niche.tiff")

