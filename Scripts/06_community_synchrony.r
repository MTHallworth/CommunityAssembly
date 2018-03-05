
A <- readRDS("Results/abund_simslist.rds")

str(A,1)

N <- A$N

AbundData <- readRDS("Data/AbunData.rds")

TempVar <- PrecipVar <- ThermVar <- combinedVar <- rep(NA,242)
for(n in 1:242){
TempVar[n] <- var(AbundData$JuneTemp[n,])
PrecipVar[n] <- var(AbundData$Precip[n,])
ThermVar[n] <- var(AbundData$ThermalSum[n,])
combinedVar[n] <- var(c(AbundData$JuneTemp[n,],AbundData$Precip[n,],AbundData$ThermalSum[n,]))
}

synch <- array(NA,c(4200,242))
for(i in 1:4200){
for(n in 1:242){
varTotN <- rep(NA,16)
SppVar <- rep(NA,13)

for(y in 1:16){
varTotN[y] <- sum(N[i,n,y,])
}
for(s in 1:13){
SppVar[s] <- sqrt(var(N[i,n,,s]))
}
synch[i,n] <- var(varTotN)/(sum(SppVar))^2
}
}


mean.sync <- apply(synch,2,mean)
SE.sync <- apply(synch,2,quantile,probs = c(0.025,0.975))
SE.sync <- SE.sync[,order(mean.sync,decreasing = TRUE)]
mean.sync <- mean.sync[order(mean.sync,decreasing = TRUE)]
#tiff("Figures/CommunityWideSync.tiff", width = 3600, height = 3600, res = 600, units = "px")
#par(bty = "l")
#plot(mean.sync,
#     ylim = c(0,1),
#     #xlim = c(1,242),
#     pch = 20,
#     ylab = "Community-wide synchrony",
#     xlab = "Survey site",
#     yaxt = "n")
#axis(2,las = 2)
#segments(x0 = 1:242,y0 = SE.sync[1,],
#         x1 = 1:242,y1 = SE.sync[2,],
#         col = "gray")
#points(1:242,mean.sync, pch = 20)
#dev.off()
system("open Figures/CommunityWideSync.tiff")


str(phi.s)
var(phi.s[,,1])
synch <- array(NA,c(4200,242))
for(i in 1:4200){
for(n in 1:242){
varTotN <- rep(NA,16)
SppVar <- rep(NA,13)

for(s in 1:13){
SppVar[s] <- var(phi.s[i,,s])
}
synch[i,n] <- var(varTotN)/(sum(SppVar))^2
}
}

mean.sync <- apply(synch,2,mean)
SE.sync <- apply(synch,2,quantile,probs = c(0.025,0.975))
SE.sync <- SE.sync[,order(mean.sync,decreasing = TRUE)]
mean.sync <- mean.sync[order(mean.sync,decreasing = TRUE)]

