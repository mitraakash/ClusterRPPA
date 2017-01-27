---
title: "Real Data for Clustering"
author: "Akash Mitra"
date: "1/18/2017"
output: R Script
---
  
  ```{r setup, include=FALSE}

rm(list=ls(all=TRUE))
library(cluster)
library(BiocInstaller)
library(biocLite("ConsensusClusterPlus"))
library(mclust)
library(flashClust)
library(FactoMineR)
library(readxl)
library(ggplot2)
library(modeltools)
library(heatmaply)
library(ConsensusClusterPlus)

####Reading data
dat <- read.csv("datafile1.csv", header=T)
d <- as.matrix.data.frame(dat[,-c(1,2,3,4,5,6,7,8)])

## Selected 80% item resampling (pItem), 80% gene resampling (pFeature), a maximum evalulated k of 6 so that cluster counts of 2,3,4,5,6 are evaluated (maxK), 50 resamplings (reps), agglomerative heirarchical clustering algorithm (clusterAlg) upon 1- Pearson correlation distances (distance), gave our output a title (title), and opted to have graphical results written to png files. We also used a specific random seed so that this example is repeatable (seed).
results = ConsensusClusterPlus(d,maxK=10,reps=100,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson", seed=100,plot="png", ml=NULL)


#For .example, the top five rows and columns of results for k=2:
results[[5]][["consensusMatrix"]][1:5,1:5]

#consensusTree - hclust object
results[[2]][["consensusTree"]]

#consensusClass - the sample classifications
results[[5]][["consensusClass"]][1:148]
results[[10]][["consensusClass"]][1:148]
plot(results[[10]][["consensusClass"]][1:148])

###Continue working from here
results[[5]][["ml"]][1:125]

#ml - consensus matrix result
#clrs - colors for cluster


####Hclust for Optimizing kmeans provided here  
icl = calcICL(results)   
icl[["clusterConsensus"]] 
icl[["itemConsensus"]][2:127,]
icl$itemConsensus$item

###PLOT - Outdated
resICl <- calcICL(results, title="Hclust for Optimizing Kmeans")

heatmap(resICl[["clusterConsensus"]])

heatmap(results[[2]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=2 cluster",  cexRow = 0.6)

heatmap(results[[3]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=3 cluster",  cexRow = 0.6)

heatmap(results[[4]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=4 cluster",  cexRow = 0.6)

heatmap(results[[5]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=5 cluster",  cexRow = 0.6)

heatmap(results[[6]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=6 cluster",  cexRow = 0.6)

heatmap(results[[7]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=7 cluster",  cexRow = 0.6)

heatmap(results[[8]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=8 cluster",  cexRow = 0.6)

heatmap(results[[9]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=9 cluster",  cexRow = 0.6)

heatmap(results[[10]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=10 cluster",  cexRow = 0.6)

###Above 5 heatmaps generated from using different k values

###Attempting to use graphics/colors



# creates a own color palette from red to green
column_annotation <- colorRampPalette(c("red", "yellow", "green"))(n = 148) #working with n=148


#column_annotation <- as.matrix(column_annotation)
heatmap(tempdat, col=palette, ColSideColors=column_annotation)


#Coloring pallette from red to green. Consensus values range from 0 (never clustered together) to 1 (always clustered together) marked by red to green.
palette <- colorRampPalette(c('#FF0000','#10FF00'))(256)
heatmap(results[[8]][["consensusMatrix"]] , Colv=F , labRow = icl$itemConsensus$item , labCol = icl$itemConsensus$item , scale='none', col = palette, cexRow = 0.4, cexCol = 0.4, ColSideColors=column_annotation)

###Best Graphic so far ^


tempdat <- results[[8]][["consensusMatrix"]]
tempdatclass <- results[[2]][["consensusClass"]]
color.map <- function(tempdatclass) { if (tempdatclass == "1") "#FF0000" else "#0000FF" }
heatcolor <- unlist(lapply(results[[2]][["consensusMatrix"]], color.map))
###Need to work with heatcolor

heatmap(tempdat, Colv=F , labRow = icl$itemConsensus$item , labCol = icl$itemConsensus$item , scale='none', col = palette, cexRow = 0.4, cexCol = 0.4, ColSideColors=heatcolor)

###Start from here
plot(ecdf(results[[2]][["consensusMatrix"]]), col = 'light blue', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[3]][["consensusMatrix"]]), add=T, col='red', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[4]][["consensusMatrix"]]), add=T, col = 'green', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[5]][["consensusMatrix"]]), add=T, col = 'dark blue', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[6]][["consensusMatrix"]]), add=T, col = 'yellow', xlab ="consensus index", ylab = "cdf")
legend(0,1, c('k=2','k=3','k=4','k=5', 'k=6') , lty=c(1,1,1,1,1), lwd=c(2.5,2.5,2.5,2.5,2.5),col=c("light blue","red","green","dark blue","yellow"))

###Creating Delta Area plot for Different K means - Continue from here on 1/18
library(deltaPlotR)




```


```{r}
rm(list=ls(all=TRUE))
library(cluster)
library(BiocInstaller)
library(biocLite("ConsensusClusterPlus"))
library(mclust)
library(flashClust)
library(FactoMineR)
library(readxl)
library(ggplot2)
library(modeltools)
library(heatmaply)

####Reading data
dat <- read.csv("datafile2.csv", header =T)
d <- as.matrix.data.frame(dat[,-c(1,2,3,4,5,6,7)])


## Selected 80% item resampling (pItem), 80% gene resampling (pFeature), a maximum evalulated k of 6 so that cluster counts of 2,3,4,5,6 are evaluated (maxK), 50 resamplings (reps), agglomerative heirarchical clustering algorithm (clusterAlg) upon 1- Pearson correlation distances (distance), gave our output a title (title), and opted to have graphical results written to png files. We also used a specific random seed so that this example is repeatable (seed).
results = ConsensusClusterPlus(d,maxK=12,reps=100,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson", seed=100,plot="png", ml=NULL)


#For .example, the top five rows and columns of results for k=2:
results[[5]][["consensusMatrix"]][1:5,1:5]

#consensusTree - hclust object
results[[2]][["consensusTree"]]

#consensusClass - the sample classifications
results[[5]][["consensusClass"]][1:178]
results[[10]][["consensusClass"]][1:178]
plot(results[[10]][["consensusClass"]][1:178])

###Continue working from here
results[[5]][["ml"]][1:125]

#ml - consensus matrix result
#clrs - colors for cluster


####Hclust for Optimizing kmeans provided here  
icl = calcICL(results)   
icl[["clusterConsensus"]] 
icl[["itemConsensus"]][2:127,]
icl$itemConsensus$item

###PLOT - Outdated
resICl <- calcICL(results, title="Hclust for Optimizing Kmeans")

heatmap(resICl[["clusterConsensus"]])

heatmap(results[[2]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=2 cluster",  cexRow = 0.6)

heatmap(results[[3]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=3 cluster",  cexRow = 0.6)

heatmap(results[[4]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=4 cluster",  cexRow = 0.6)

heatmap(results[[5]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=5 cluster",  cexRow = 0.6)

heatmap(results[[6]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=6 cluster",  cexRow = 0.6)

heatmap(results[[7]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=7 cluster",  cexRow = 0.6)

heatmap(results[[8]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=8 cluster",  cexRow = 0.6)

heatmap(results[[9]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=9 cluster",  cexRow = 0.6)

heatmap(results[[10]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=10 cluster",  cexRow = 0.6)

###Above 5 heatmaps generated from using different k values

###Attempting to use graphics/colors



# creates a own color palette from red to green
column_annotation <- colorRampPalette(c("red", "yellow", "green"))(n = 179) #working with n=178


#column_annotation <- as.matrix(column_annotation)
heatmap(tempdat, col=palette, ColSideColors=column_annotation)


#Coloring pallette from red to green. Consensus values range from 0 (never clustered together) to 1 (always clustered together) marked by red to green.
palette <- colorRampPalette(c('#FF0000','#10FF00'))(256)
heatmap(results[[10]][["consensusMatrix"]] , Colv=F , labRow = icl$itemConsensus$item , labCol = icl$itemConsensus$item , scale='none', col = palette, cexRow = 0.4, cexCol = 0.4, ColSideColors=column_annotation)

###Best Graphic so far ^


tempdat <- results[[10]][["consensusMatrix"]]
tempdatclass <- results[[2]][["consensusClass"]]
color.map <- function(tempdatclass) { if (tempdatclass == "1") "#FF0000" else "#0000FF" }
heatcolor <- unlist(lapply(results[[2]][["consensusMatrix"]], color.map))
###Need to work with heatcolor

heatmap(tempdat, Colv=F , labRow = icl$itemConsensus$item , labCol = icl$itemConsensus$item , scale='none', col = palette, cexRow = 0.4, cexCol = 0.4, ColSideColors=heatcolor)

###Start from here
plot(ecdf(results[[2]][["consensusMatrix"]]), col = 'light blue', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[3]][["consensusMatrix"]]), add=T, col='red', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[4]][["consensusMatrix"]]), add=T, col = 'green', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[5]][["consensusMatrix"]]), add=T, col = 'dark blue', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[6]][["consensusMatrix"]]), add=T, col = 'yellow', xlab ="consensus index", ylab = "cdf")
legend(0,1, c('k=2','k=3','k=4','k=5', 'k=6') , lty=c(1,1,1,1,1), lwd=c(2.5,2.5,2.5,2.5,2.5),col=c("light blue","red","green","dark blue","yellow"))

###Creating Delta Area plot for Different K means - Continue from here on 1/18
library(deltaPlotR)




```


```{r}
rm(list=ls(all=TRUE))
library(cluster)
library(BiocInstaller)
library(biocLite("ConsensusClusterPlus"))
library(mclust)
library(flashClust)
library(FactoMineR)
library(readxl)
library(ggplot2)
library(modeltools)
library(heatmaply)

####Reading data
dat <- read.csv("datafile3.csv", header =T)
d <- as.matrix.data.frame(dat[,-c(1,2,3,4,5,6,7)])


## Selected 80% item resampling (pItem), 80% gene resampling (pFeature), a maximum evalulated k of 6 so that cluster counts of 2,3,4,5,6 are evaluated (maxK), 50 resamplings (reps), agglomerative heirarchical clustering algorithm (clusterAlg) upon 1- Pearson correlation distances (distance), gave our output a title (title), and opted to have graphical results written to png files. We also used a specific random seed so that this example is repeatable (seed).
results = ConsensusClusterPlus(d,maxK=10,reps=100,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson", seed=100,plot="png", ml=NULL)


#For .example, the top five rows and columns of results for k=2:
results[[5]][["consensusMatrix"]][1:5,1:5]

#consensusTree - hclust object
results[[2]][["consensusTree"]]

#consensusClass - the sample classifications
results[[5]][["consensusClass"]][1:175]
results[[10]][["consensusClass"]][1:175]
plot(results[[10]][["consensusClass"]][1:175])

###Continue working from here
results[[5]][["ml"]][1:125]

#ml - consensus matrix result
#clrs - colors for cluster


####Hclust for Optimizing kmeans provided here  
icl = calcICL(results)   
icl[["clusterConsensus"]] 
icl[["itemConsensus"]][2:127,]
icl$itemConsensus$item

###PLOT - Outdated
resICl <- calcICL(results, title="Hclust for Optimizing Kmeans")

heatmap(resICl[["clusterConsensus"]])

heatmap(results[[2]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=2 cluster",  cexRow = 0.6)

heatmap(results[[3]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=3 cluster",  cexRow = 0.6)

heatmap(results[[4]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=4 cluster",  cexRow = 0.6)

heatmap(results[[5]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=5 cluster",  cexRow = 0.6)

heatmap(results[[6]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=6 cluster",  cexRow = 0.6)

heatmap(results[[7]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=7 cluster",  cexRow = 0.6)

heatmap(results[[8]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=8 cluster",  cexRow = 0.6)

heatmap(results[[9]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=9 cluster",  cexRow = 0.6)

heatmap(results[[10]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=10 cluster",  cexRow = 0.6)

###Above 5 heatmaps generated from using different k values

###Attempting to use graphics/colors

tempdat <- results[[8]][["consensusMatrix"]]

# creates a own color palette from red to green
column_annotation <- colorRampPalette(c("red", "yellow", "green"))(n = 176) #working with n=176


#column_annotation <- as.matrix(column_annotation)
heatmap(tempdat, col=palette, ColSideColors=column_annotation)


#Coloring pallette from red to green. Consensus values range from 0 (never clustered together) to 1 (always clustered together) marked by red to green.
palette <- colorRampPalette(c('#FF0000','#10FF00'))(256)
heatmap(results[[8]][["consensusMatrix"]] , Colv=F , labRow = icl$itemConsensus$item , labCol = icl$itemConsensus$item , scale='none', col = palette, cexRow = 0.4, cexCol = 0.4, ColSideColors=column_annotation)

###Best Graphic so far ^



tempdatclass <- results[[2]][["consensusClass"]]
color.map <- function(tempdatclass) { if (tempdatclass == "1") "#FF0000" else "#0000FF" }
heatcolor <- unlist(lapply(results[[2]][["consensusMatrix"]], color.map))
###Need to work with heatcolor

heatmap(tempdat, Colv=F , labRow = icl$itemConsensus$item , labCol = icl$itemConsensus$item , scale='none', col = palette, cexRow = 0.4, cexCol = 0.4, ColSideColors=heatcolor)

###Start from here
plot(ecdf(results[[2]][["consensusMatrix"]]), col = 'light blue', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[3]][["consensusMatrix"]]), add=T, col='red', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[4]][["consensusMatrix"]]), add=T, col = 'green', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[5]][["consensusMatrix"]]), add=T, col = 'dark blue', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[6]][["consensusMatrix"]]), add=T, col = 'yellow', xlab ="consensus index", ylab = "cdf")
legend(0,1, c('k=2','k=3','k=4','k=5', 'k=6') , lty=c(1,1,1,1,1), lwd=c(2.5,2.5,2.5,2.5,2.5),col=c("light blue","red","green","dark blue","yellow"))

###Creating Delta Area plot for Different K means - Continue from here on 1/18
library(deltaPlotR)




```


```{r}
rm(list=ls(all=TRUE))
library(cluster)
library(BiocInstaller)
library(biocLite("ConsensusClusterPlus"))
library(mclust)
library(flashClust)
library(FactoMineR)
library(readxl)
library(ggplot2)
library(modeltools)
library(heatmaply)

####Reading data
dat <- read.csv("datafile4.csv", header =T)
d <- as.matrix.data.frame(dat[,-c(1,2,3,4,5,6,7)])


## Selected 80% item resampling (pItem), 80% gene resampling (pFeature), a maximum evalulated k of 6 so that cluster counts of 2,3,4,5,6 are evaluated (maxK), 50 resamplings (reps), agglomerative heirarchical clustering algorithm (clusterAlg) upon 1- Pearson correlation distances (distance), gave our output a title (title), and opted to have graphical results written to png files. We also used a specific random seed so that this example is repeatable (seed).
results = ConsensusClusterPlus(d,maxK=10,reps=100,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson", seed=100,plot="png", ml=NULL)


#For .example, the top five rows and columns of results for k=2:
results[[5]][["consensusMatrix"]][1:5,1:5]

#consensusTree - hclust object
results[[2]][["consensusTree"]]

#consensusClass - the sample classifications
results[[5]][["consensusClass"]][1:148]
results[[10]][["consensusClass"]][1:148]
plot(results[[10]][["consensusClass"]][1:148])

###Continue working from here
results[[5]][["ml"]][1:125]

#ml - consensus matrix result
#clrs - colors for cluster


####Hclust for Optimizing kmeans provided here  
icl = calcICL(results)   
icl[["clusterConsensus"]] 
icl[["itemConsensus"]][2:127,]
icl$itemConsensus$item

###PLOT - Outdated
resICl <- calcICL(results, title="Hclust for Optimizing Kmeans")

heatmap(resICl[["clusterConsensus"]])

heatmap(results[[2]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=2 cluster",  cexRow = 0.6)

heatmap(results[[3]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=3 cluster",  cexRow = 0.6)

heatmap(results[[4]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=4 cluster",  cexRow = 0.6)

heatmap(results[[5]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=5 cluster",  cexRow = 0.6)

heatmap(results[[6]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=6 cluster",  cexRow = 0.6)

heatmap(results[[7]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=7 cluster",  cexRow = 0.6)

heatmap(results[[8]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=8 cluster",  cexRow = 0.6)

heatmap(results[[9]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=9 cluster",  cexRow = 0.6)

heatmap(results[[10]][["consensusMatrix"]], scale = c("row", "column", "none"), na.rm = TRUE, margins = c(3, 3) , main = "k=10 cluster",  cexRow = 0.6)

###Above 5 heatmaps generated from using different k values

###Attempting to use graphics/colors



# creates a own color palette from red to green
column_annotation <- colorRampPalette(c("red", "yellow", "green"))(n = 148) #working with n=148


#column_annotation <- as.matrix(column_annotation)
heatmap(tempdat, col=palette, ColSideColors=column_annotation)


#Coloring pallette from red to green. Consensus values range from 0 (never clustered together) to 1 (always clustered together) marked by red to green.
palette <- colorRampPalette(c('#FF0000','#10FF00'))(256)
heatmap(results[[8]][["consensusMatrix"]] , Colv=F , labRow = icl$itemConsensus$item , labCol = icl$itemConsensus$item , scale='none', col = palette, cexRow = 0.4, cexCol = 0.4, ColSideColors=column_annotation)

###Best Graphic so far ^


tempdat <- results[[2]][["consensusMatrix"]]
tempdatclass <- results[[2]][["consensusClass"]]
color.map <- function(tempdatclass) { if (tempdatclass == "1") "#FF0000" else "#0000FF" }
heatcolor <- unlist(lapply(results[[2]][["consensusMatrix"]], color.map))
###Need to work with heatcolor

heatmap(tempdat, Colv=F , labRow = icl$itemConsensus$item , labCol = icl$itemConsensus$item , scale='none', col = palette, cexRow = 0.4, cexCol = 0.4, ColSideColors=heatcolor)

###Start from here
plot(ecdf(results[[2]][["consensusMatrix"]]), col = 'light blue', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[3]][["consensusMatrix"]]), add=T, col='red', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[4]][["consensusMatrix"]]), add=T, col = 'green', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[5]][["consensusMatrix"]]), add=T, col = 'dark blue', xlab ="consensus index", ylab = "cdf")
plot(ecdf(results[[6]][["consensusMatrix"]]), add=T, col = 'yellow', xlab ="consensus index", ylab = "cdf")
legend(0,1, c('k=2','k=3','k=4','k=5', 'k=6') , lty=c(1,1,1,1,1), lwd=c(2.5,2.5,2.5,2.5,2.5),col=c("light blue","red","green","dark blue","yellow"))

###Creating Delta Area plot for Different K means - Continue from here on 1/18
library(deltaPlotR)





```
