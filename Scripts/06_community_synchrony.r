
A <- readRDS("DynAbun_simslist_2018-03-16.rds")

nsims <- dim(A$N)[[1]]

synch <- array(NA,c(nsims,242))

for(i in 1:nsims){
for(n in 1:242){
varTotN <- rep(NA,16)
SppVar <- rep(NA,13)

for(y in 1:16){
varTotN[y] <- sum(A$N[i,n,y,])
}
for(s in 1:13){
SppVar[s] <- sqrt(var(A$N[i,n,,s]))
}
synch[i,n] <- var(varTotN)/(sum(SppVar))^2
}
}

saveRDS(synch,"synch.rds")

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


######################################################################################################
#
#
# Species abundance distribution using the GAMBIN model
#
#
######################################################################################################
library(gambin, lib.loc = "../CostOfRepro/Packages")

Alphas <- vector("list",16)
for(y in 1:16){
cat(paste0("Currently analyzing year ",y,"\n"))
alpha <- array(NA,c(1000,242))
cat(paste0("Started at ", Sys.time(),"\n"))
a <- Sys.time()
for(i in 1:242){
TEMP <- A$N[,i,y,]
TEMP <- TEMP[rowSums(TEMP)>15,]
if(nrow(TEMP)>1000){
TEMP <- TEMP[sample(1000),]
}
else{TEMP <- TEMP[sample(1:nrow(TEMP),1000,replace = TRUE),]
}
alpha[,i] <- apply(TEMP, 1, FUN = function(x){fit_abundances(x,subsample = 15)$alpha})
}
Alphas[[y]]<-alpha
cat(paste0("Finished - it took ",Sys.time()-a," mins\n"))
}

saveRDS(Alphas,"Alphas_dynabun.rds")



