library(caret)
library(clusteval) 
library(cluster)
library(corrplot)
library(d3heatmap)
library(fpc)
library(gplots)
library(NMF)

#--------------------------------
# Clustering Functions
#--------------------------------

# colors a htclus (hc) by "K" a number of clusters
dendColorOnK <- function(hc, kclust) {
  labelColors = c("#ff3333","#33ccff","#33ff33","#cc33ff","#ffff33","#ff9933")
  #labelColors = palette(rainbow(kclust))
  clusMember = cutree(hc, kclust)
  colClus <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
      #print(clusMember[which(names(clusMember) == a$label)]) 
    }
    n
  }
  hcd = as.dendrogram(hc)
  clusDendro = dendrapply(hcd, colClus)
  return(clusDendro)
}

# colors according to cluster
dendColorOnClus <- function(hc, sampInfo, nameCol, clusCol) {
  labelColors = matrix(
    c(1,2,3,4,5,6,"#ff3333","#33ccff","#33ff33","#cc33ff","#ffff33","#ff9933"), nrow=6, ncol=2
  )
  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      # This removes the .fasta and then looks up the column in "t"
      haystack = gsub("\\..*","",a$label)
      #haystack = a$label
      indicies <- haystack == sampInfo[nameCol]
      clusNum = sampInfo[clusCol][indicies]
      labCol <- labelColors[as.integer(clusNum),2]
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
  }
  hcd = as.dendrogram(hc)
  clusDendro = dendrapply(hcd, colLab)
  return(clusDendro)
}

# goes from samp 1. samp 2. distance to matrix
nameDistToMatrix <- function(x){
  x.names <- sort(unique(c(x[[1]], x[[2]])))
  x.dist <- matrix(0, length(x.names), length(x.names))
  dimnames(x.dist) <- list(x.names, x.names)
  x.ind <- rbind(cbind(match(x[[1]], x.names), match(x[[2]], x.names)), cbind(match(x[[2]], x.names), match(x[[1]], x.names)))
  x.dist[x.ind] <- rep(x[[3]], 2)
  return(x.dist)
}

# Jaccard simulairity between two clusters
clusterSimularityJaccard <- function(clusteringsTable){
  simularityList = matrix(ncol=3)
  for( i in 2:ncol(clusteringsTable) ){
    for( n in i:ncol(clusteringsTable) ){
      colnamei <- colnames(clusteringsTable[i])
      colnamen <- colnames(clusteringsTable[n])
      jacdist <- as.double(clusteval::cluster_similarity(clusteringsTable[[i]], clusteringsTable[[n]], similarity = c("jaccard"), method = "independence"))
      nameDistRow <- cbind(colnamei,colnamen,jacdist)
      simularityList <- rbind(simularityList, nameDistRow)
    }
  }
  return(simularityList)
}

# creates a table of simulariteis 
clusterSimularityTable <- function(clusteringsTable){
  simularityList = matrix(ncol=6)
  for( i in 2:ncol(clusteringsTable) ){
    for( n in i:ncol(clusteringsTable) ){
      colnamea <- colnames(clusteringsTable[i])
      colnameb <- colnames(clusteringsTable[n])
      jac <- cluster_similarity(clusteringsTable[[i]], clusteringsTable[[n]], similarity = c("jaccard"), method = "independence")
      rand <- cluster_similarity(clusteringsTable[[i]], clusteringsTable[[n]], similarity = c("rand"), method = "independence")
      wilcox <- wilcox.test(clusteringsTable[,i],clusteringsTable[,n])
      wilcox_p_value <- wilcox$p.value
      wilcox_statistic <- wilcox$statistic
      nameDistRow <- cbind(colnamea,colnameb,jac,rand,wilcox_p_value,wilcox_statistic )
      simularityList <- rbind(simularityList, nameDistRow)
    }
  }
  return(simularityList)
}


cool1=c("#cc33ff","#33ccff", "#33ff33")
cool2=c("#cc33ff","#33ccff", "#33ff33","#00FFFF")

col<- colorRampPalette(c("white", "#33ccff"))(100)
cola<- colorRampPalette(c("#33ccff", "white"))(100)
colb<- colorRampPalette(c("#33ccff","white","#ff3333"))(200)
col1<- colorRampPalette(c("#ff3333","white","#33ff33"))(200)
col2<- c("#33ff33", "white")
col3<- c("#ff3333")
col4<- c("#ff9933")
col5<- c("#cc33ff","#33ff33","#33ccff","#ff3333")

r <- read.csv(“otu_table.csv", as.is=TRUE)
x <- read.table("Atacama_mash_dists.txt", as.is=TRUE)
#x1 <- read.table("Untrim_Atacama_Dists.txt", as.is=TRUE)
sra <- read.csv(“metadata.csv", as.is=TRUE)
f1 <- read.csv(“clusters.K3.csv",  as.is=TRUE)
f2 <- read.csv("clusters.K4.csv",  as.is=TRUE)
f1j <- read.csv("k3_jaccard_sim_between_groupings.csv", header = TRUE,  as.is=TRUE)
f2j <- read.csv("k4_jaccard_sim_between_groupings.csv", header = TRUE,  as.is=TRUE)

# swtich jaccard matrix of clusters to matrix format
f1j.dist <- nameDistToMatrix(f1j)
f2j.dist <- nameDistToMatrix(f2j)

# generate a cluster simularity table for k=3 and k=4
f1.meas <- clusterSimularityTable(f1)
f2.meas <- clusterSimularityTable(f2)
###write.csv(f1.meas, file = "k3_cluster_Simularity_Table.csv")
###write.csv(f2.meas, file = "k4_cluster_Simularity_Table.csv")

# generate a cluster simularity martix via jaccard for k=3 and k=4
f1.jac <- clusterSimularityJaccard(f1)
f2.jac <- clusterSimularityJaccard(f2)
###write.csv(f1.jac, file = "k3_jaccard_sim_between_groupings.csv")
###write.csv(f2.jac, file = "k4_jaccard_sim_between_groupings.csv")

# import and transform mash distances
x.dist <- nameDistToMatrix(x)
###hc = hclust(dist(x.dist), method="mcquitty")
###cutdendo = dendColorOnK(hc, 4)
AtacamaSamples = x.dist
dianaObj <- diana(AtacamaSamples)

#Begin Silohette
comp.sra <- sra[,1:7]
comp.sra$loc <- sra[,27]

trans = preProcess(comp.sra[,3:8],
                   method=c("BoxCox", "center", 
                            "scale", "pca"))

PC = predict(trans, comp.sra[,3:8])
k3PCAkmeans <- kmeans(PC,3)
k4PCAkmeans <- kmeans(PC,4)
k3kmeans <- kmeans(as.dist(x.dist),3)
k4kmeans <- kmeans(as.dist(x.dist),4)

OTUdata <- log(r[,2:61])
#OTUdata <- r[,2:61]
OTUdata.clean <- do.call(data.frame,lapply(OTUdata, function(x) replace(x, is.infinite(x),0)))
OTUdata.clean <- OTUdata.clean[rowSums(OTUdata.clean) > 0,]
OTUdata.k3kmeans <- kmeans(OTUdata.clean[2:60],3)
OTUdata.k4kmeans <- kmeans(OTUdata.clean[2:60],4)

OTU.trans = preProcess(OTUdata.clean[2:60],
                               method=c("BoxCox", "center", 
                                        "scale", "pca"))
OTUdata.PC = predict(OTU.trans, OTUdata.clean[2:60])

#For number of clusters
#wss <- (nrow(OTUdata.clean[2:60])-1)*sum(apply(OTUdata.clean[2:60],2,var))
#for (i in 2:15) wss[i] <- sum(kmeans(OTUdata.clean[2:60],
#                                     centers=i)$withinss)
#plot(1:15, wss, type="b", xlab="Number of Clusters",
#     ylab="Within groups sum of squares")

OTUdata.k3kmeans <- kmeans(OTUdata.PC,3)
OTUdata.k4kmeans <- kmeans(OTUdata.PC,4)

OTUdata.k3nmf <- nmf(OTUdata.clean[2:60],3,'nsNMF')
OTUdata.k4nmf <- nmf(OTUdata.clean[2:60],4,'nsNMF')

dianak3 <- cutree(as.hclust(dianaObj),3)
dianak4 <- cutree(as.hclust(dianaObj),4)

mcquittyk3 <- cutree(hclust(as.dist(x.dist),method="mcquitty"),3)
mcquittyk4 <- cutree(hclust(as.dist(x.dist),method="mcquitty"),4)

#png(filename = "silhouettes.png", width = 1024, height = 1792)
#par( mfrow = c( 7, 2 ) )

OTUdata.k3kmeans.dissE <- daisy(OTUdata.PC)
OTUdata.k3kmeans.dE2   <- OTUdata.k3kmeans.dissE^2
OTUdata.k3kmeans.sk2   <- silhouette(OTUdata.k3kmeans$cluster, OTUdata.k3kmeans.dE2)
plot(OTUdata.k3kmeans.sk2, col=cool1, main="KMEANS PCA OTU Abundance, K=3", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

OTUdata.k4kmeans.dissE <- daisy(OTUdata.PC)
OTUdata.k4kmeans.dE2   <- OTUdata.k4kmeans.dissE^2
OTUdata.k4kmeans.sk2   <- silhouette(OTUdata.k4kmeans$cluster, OTUdata.k4kmeans.dE2)
plot(OTUdata.k4kmeans.sk2, col=cool2, main="KMEANS PCA OTU Abundance, K=4", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

k3PCAkmeans.dissE <- daisy(PC)
k3PCAkmeans.dE2   <- k3PCAkmeans.dissE^2
k3PCAkmeans.sk2   <- silhouette(k3PCAkmeans$cluster, k3PCAkmeans.dE2)
plot(k3PCAkmeans.sk2, col=cool1, main="KMEANS PCA ENV, K=3", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

k4PCAkmeans.dissE <- daisy(PC)
k4PCAkmeans.dE2   <- k4PCAkmeans.dissE^2
k4PCAkmeans.sk2   <- silhouette(k4PCAkmeans$cluster, k4PCAkmeans.dE2)
plot(k4PCAkmeans.sk2, col=cool2, main="KMEANS PCA ENV, K=4", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

k3kmeans.dissE <- as.dist(x.dist)
k3kmeans.dE2   <- k3kmeans.dissE^2
k3kmeans.sk2   <- silhouette(k3kmeans$cl, k3kmeans.dE2)
plot(k3kmeans.sk2, col=cool1, main="KMEANS DIST, K=3", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

k4kmeans.dissE <- as.dist(x.dist)
k4kmeans.dE2   <- k4kmeans.dissE^2
k4kmeans.sk2   <- silhouette(k4kmeans$cl, k4kmeans.dE2)
plot(k4kmeans.sk2, col=cool2, main="KMEANS DIST, K=4", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

dianaObj.k3.sil   <- silhouette(dianak3, dianaObj$diss)
plot(dianaObj.k3.sil, col=cool1, main="DIANA CUTREE, K=3", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

dianaObj.k4.sil   <- silhouette(dianak4, dianaObj$diss)
plot(dianaObj.k4.sil, col=cool2, main="DIANA CUTREE, K=4", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

mcquittyk3.dissE <- as.dist(x.dist)
mcquittyk3.dE2 <- mcquittyk3.dissE^2
mcquittyk3.sil <- silhouette(mcquittyk3, mcquittyk3.dE2)
plot(mcquittyk3.sil, col=cool1, main="MCQUITTY CUTREE, K=3", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUEplo)

mcquittyk4.dissE <- as.dist(x.dist)
mcquittyk4.dE2 <- mcquittyk4.dissE^2
mcquittyk4.sil <- silhouette(mcquittyk4, mcquittyk4.dE2)
plot(mcquittyk4.sil, col=cool2, main="MCQUITTY CUTREE, K=3", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

OTUdata.k3kmeans.dissE <- daisy(OTUdata.clean)
OTUdata.k3kmeans.dE2   <- OTUdata.k3kmeans.dissE^2
OTUdata.k3kmeans.sk2   <- silhouette(OTUdata.k3kmeans$cl, OTUdata.k3kmeans.dE2)
plot(OTUdata.k3kmeans.sk2, col=cool1, main="OTU KMEANS, K=3", cex.axis=2, cex.names=2, cex.lab=0.01, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

OTUdata.k4kmeans.dissE <- daisy(OTUdata.clean)
OTUdata.k4kmeans.dE2   <- OTUdata.k4kmeans.dissE^2
OTUdata.k4kmeans.sk2   <- silhouette(OTUdata.k4kmeans$cl, OTUdata.k4kmeans.dE2)
plot(OTUdata.k4kmeans.sk2, col=cool2, main="OTU KMEANS, K=4", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

OTUdata.k3nmf.sil <- silhouette(OTUdata.k3nmf)
plot(OTUdata.k3nmf.sil, col=cool1, main="NMF OTU LOG, K=3", cex.axis=2, cex.names=2, cex.lab=0.01,  do.n.k = FALSE, do.clus.stat = TRUE)

OTUdata.k4nmf.sil <- silhouette(OTUdata.k4nmf)
plot(OTUdata.k4nmf.sil, col=cool2, main="NMF OTU LOG, K=4", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)
dev.off()

#!png(filename = "environmentalPCA.png", width = 256, height = 512)
par( mfrow = c( 2, 1 ) )
plotcluster(PC, k3PCAkmeans$cluster, col=col5[k3PCAkmeans$cluster], pch=16, method="dc")
plotcluster(PC, k4PCAkmeans$cluster, col=col5[k4PCAkmeans$cluster], pch=16, method="dc")
#!dev.off()

with(comp.sra, pairs(comp.sra[,3:8], col=col5[k3PCAkmeans$cluster])) 
with(comp.sra, pairs(comp.sra[,3:8], col=col5[k4PCAkmeans$cluster])) 

plot(dianaObj, which.plots=1, col=col2)

# Corr plot for k=3 and k=4 
corrplot(f1j.dist, order="hclust", method="color", col=col1, cl.lim=c(0,1), is.corr=FALSE, hclust.method="centroid", tl.col="black", addgrid.col="white", addCoef.col="black", type="upper")
corrplot(f2j.dist, order="hclust", method="color", col=col1, cl.lim=c(0,1), is.corr=FALSE, hclust.method="centroid", tl.col="black", addgrid.col="white", addCoef.col="black", type="upper")

# Heatmap clustering functions
AtacamaSamples.funa <- function(x) hclust(dist(x), method="mcquitty")
AtacamaSamples.funb <- function(x) diana(dist(x))

#heatmaps
heatmap(AtacamaSamples, col = cola, symm = TRUE, hclustfun = AtacamaSamples.funa, main = "All distances from Sample to Sample, clustered by mcquitty")
heatmap(AtacamaSamples, col = cola, symm = TRUE, hclustfun = AtacamaSamples.funb, main = "All distances from Sample to Sample, clustered by diana")


# Location colored dendograms
par( mfrow = c( 2, 1 ) )
cutdendo = dendColorOnClus(dianaObj,sra,1,27)
plot(cutdendo, ylab="Distance", main = "Sampling Location by Site, K=6")
cutdendo = dendColorOnClus(dianaObj,f1,1,8)
plot(cutdendo, ylab="Distance", main = "Sampling Location by Cardinal Direction, K=3")

# Dendograms colored according to group, k=3
cutdendo = dendColorOnClus(dianaObj,f1,1,2)
plot(cutdendo, ylab="Distance", main = "NMF OTU, K=3")
cutdendo = dendColorOnK(as.hclust(dianaObj), 3)
plot(cutdendo, ylab="Distance", main = "DIANA CUTREE, K=3")
cutdendo = dendColorOnClus(dianaObj,f1,1,4)
plot(cutdendo, ylab="Distance", main = "PH RANGE, K=3")
cutdendo = dendColorOnClus(dianaObj,f1,1,5)
plot(cutdendo, ylab="Distance", main = "NMF Subtypes with Log transform, K=3")
cutdendo = dendColorOnClus(dianaObj,f1,1,6)
plot(cutdendo, ylab="Distance", main = "RANDOM, K=3")
cutdendo = dendColorOnClus(dianaObj,f1,1,7)
plot(cutdendo, ylab="Distance", main = "RANDOM 2, K=3")
cutdendo = dendColorOnClus(dianaObj,f1,1,8)
plot(cutdendo, ylab="Distance", main = "CARDINAL DIRECTION, K=3")
cutdendo = dendColorOnClus(dianaObj,f1,1,9)
plot(cutdendo, ylab="Distance", main = "MCQUITTY CUTREE, K=3")
cutdendo = dendColorOnClus(dianaObj,f1,1,10)
plot(cutdendo, ylab="Distance", main = "KMEANS PCA ENV, K=3")
cutdendo = dendColorOnClus(dianaObj,f1,1,11)
plot(cutdendo, ylab="Distance", main = "KMEANS DIST, K=3")

# Dendograms colored according to group, k=4
cutdendo = dendColorOnClus(dianaObj,f2,1,2)
plot(cutdendo, ylab="Distance", main = "RANDOM, K=4")
cutdendo = dendColorOnClus(dianaObj,f2,1,3)
plot(cutdendo, ylab="Distance", main = "RANDOM 2, K=4")
cutdendo = dendColorOnClus(dianaObj,f2,1,4)
plot(cutdendo, ylab="Distance", main = "PH RANGE PH7 DEPTH MIX, K=4")
cutdendo = dendColorOnClus(dianaObj,f2,1,5)
plot(cutdendo, ylab="Distance", main = "NMF LOG OTU, K=4")
cutdendo = dendColorOnK(as.hclust(dianaObj), 4)
plot(cutdendo, ylab="Distance", main = "DIANA CUTREE, K=4")
cutdendo = dendColorOnClus(dianaObj,f2,1,7)
plot(cutdendo, ylab="Distance", main = "PH RANGE, K=4")
cutdendo = dendColorOnClus(dianaObj,f2,1,8)
plot(cutdendo, ylab="Distance", main = "MCQUITTY CUTREE, K=4")
cutdendo = dendColorOnClus(dianaObj,f2,1,9)
plot(cutdendo, ylab="Distance", main = "PKMEANS PCA ENV, K=4")
cutdendo = dendColorOnClus(dianaObj,f2,1,10)
plot(cutdendo, ylab="Distance", main = "KMEANS DIST, K=4")

# Anova section
# Clean up Env. Var.
sram <- log(sra[3:27])
sram <- do.call(data.frame,lapply(sram, function(x) replace(x, is.infinite(x),NaN)))
sram[,5][is.na(sram[,5])] <- 0


man.lm <- lm(sram$NMF_OTU_K3~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$DIANA_CUTREE_K3~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$NMF_LOG_OTU_K3~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$MCQUITTY_CUTREE_K3~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$KMEANS_DIST_K3~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$KMEANS_PCA_ENV_K3~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)


man.lm <- lm(sram$PH_RANGE_K3~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$RANDOM_K3~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$RANDOM_2_K3~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$CARDINAL_DIRECTION_K3~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$RANDOM_K4~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$RANDOM_2_K4~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$KMEANS_PCA_ENV_K4~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$PH_RANGE_PH7_DEPTH_MIX_K4~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$PH_RANGE_K4~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$NMF_LOG_OTU_K4~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$DIANA_CUTREE_K4~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$MCQUITTY_CUTREE_K4~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)

man.lm <- lm(sram$KMEANS_DIST_K4~sram$pH*sram$Elevation*sram$Conductivity*sram$Air_Relative_Humidity*sram$Depth, data=sram)
man.res <- aov(man.lm)
man.ll <- logLik(man.res)
summary(man.ll)
print(man.ll)
summary(man.res)
#!plot(man.res, col=col3)

dev.off()

#Heatmaps for env var relation
tab.sra <- sra[,1:7]
tab.sra$loc <- sra[,27]
tab.sra$Sample_Location_Name <- NULL

tab = tab.sra[,-1]
tab = tab
tab = data.frame(tab)
tab = do.call(data.frame,lapply(tab, function(x) replace(x, is.infinite(x),0)))
rownames(tab) <- tab.sra[,1]

tab <- tab[!rownames(tab) %in% c("SRR901092"), ]

d3heatmap(tab,Colv = NULL, Rowv = as.dendrogram(dianaObj), scale="column", col=colb, cexCol = .65, cexRow = .5)

DianaDendrogram <- as.dendrogram(dianaObj)

k3tab = f1[,-1]
k3tab = data.frame(k3tab)
k3tab.hclust <- hclust(dist(1-f1j.dist), method="centroid")
k3tab.col <- as.dendrogram(k3tab.hclust)
rownames(k3tab) <- f1[,1]
k3tab <- k3tab[!rownames(k3tab) %in% c("SRR901092"), ]


heatmap.2(as.matrix(k3tab), trace="none", key="FALSE", dendrogram="row", margins=c(10,50),
        cexCol = .8, #cexRow = .5,
        Rowv = as.dendrogram(dianaObj), Colv = k3tab.col, 
        col=c("#ffff33","#33ccff","#33ff33","#ff3333"))

k4tab = f2[,-1]
k4tab = data.frame(k4tab)
k4tab.hclust <- hclust(dist(1-f2j.dist), method="centroid")
k4tab.col <- as.dendrogram(k4tab.hclust)
rownames(k4tab) <- f2[,1]
k4tab <- k4tab[!rownames(k4tab) %in% c("SRR901092"), ]

heatmap.2(as.matrix(k4tab), trace="none", key="FALSE", dendrogram="row", margins=c(10,50),
        cexCol = .6, #cexRow = .5,
        Rowv = as.dendrogram(dianaObj), Colv = k4tab.col, 
        col=c("#ffff33","#33ccff","#33ff33","#ff3333","#cc33ff"))

