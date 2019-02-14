library(ecodist)
library(ggplot2)
library(vegan)
library(reshape2)
library(zoo)
library(readr)
update.packages(ask = FALSE)

## Chapter 2
## mds
## Data preparation

abund_data <-  read.csv(file = "PB_data_matrix.csv", row.names = 'Replicate', sep = ",", header = TRUE) 
Site  <- abund_data$Replicate
meta_table<-read.csv(file = "PB_diat_env.csv",row.names= "Replicate",sep = ",", header = TRUE)
# alldata2 <- subset(alldata2,select=-c(Site))

all_bray <- bcdist(abund_data) # Create bray - curtis dissimilarity matrix

#fauna <-  read.csv(file = "copy of Alldata_fauna.csv", row.names = 'Replicate', sep = ",", header = TRUE) # All dry biomass of species for all sites and replicates for fauna only (seperate file)
#Site  <- fauna$Site
#fauna <- subset(fauna,select=-c(Site))

#flora <-  read.csv(file = "copy of Alldata_flora.csv", row.names = 'Replicate', sep = ",", header = TRUE) # All dry biomass of species for all sites and replicates for flora only (seperate file)
#Site  <- flora$Site
#flora <- subset(flora,select=-c(Site))

#fauna_matrix <- bcdist(fauna) # Create bray - curtis dissimilarity matrix fauna
#flora_matrix <- bcdist(flora) # Create bray - curtis dissimilarity matrix flora

benthic.mds <- metaMDS(abund_data)
benthic.mds
?metaMDS
####################



#######


data.scores <- as.data.frame(scores(benthic.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(benthic.mds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data


## Plotting mds: Method 1
#
plot(benthic.mds, type = "n")
points(benthic.mds, display = "sites", cex = 0.8, pch = c(fill = Site), col = "black")
ordihull(benthic.mds, meta_table$Host_spp, col=1:3)

#
plot(benthic.mds, type = "n")
points(benthic.mds, display = "sites", cex = 0.8, pch = 21, col = c(fill = Site))
ordiellipse(benthic.mds, meta_table$Host_spp, col=1:3, draw="polygon")
ordispider(benthic.mds, meta_table$Host_spp, Site, col=1:3, label = TRUE)

# My nMDS plots
# $Host_spp
pdf("NMDSplot1.pdf")
#plot(benthic.mds, type = "n")
plot(benthic.mds$points, col = meta_table$Host_spp)
#points(benthic.mds, display = "sites", cex = 0.8, pch = c(fill = Site), col = "black")
ordiellipse(benthic.mds, meta_table$Host_spp, col=1:3, draw="polygon")
ordihull(benthic.mds, meta_table$Host_spp,col=1:3)
ordispider(benthic.mds, meta_table$Host_spp, col= 1:3, label = TRUE)
dev.off()
# $Host_size
pdf("NMDSplot.pdf")
#plot(benthic.mds$points, col = meta_table$Host_size)
plot(benthic.mds, type = "n")
#points(benthic.mds, display = "sites", cex = 0.8, pch = c(fill = Site), col = "black")
ordiellipse(benthic.mds, meta_table$Host_size, col=1:3, draw="polygon")
ordihull(benthic.mds, meta_table$Host_size,col=1:3)
ordispider(benthic.mds, meta_table$Host_size, col= 1:3, label = TRUE)
dev.off()
#
plot(benthic.mds, type = "n")
points(benthic.mds, display = "sites", cex = 0.8, pch = c(fill = Site), col = c(fill = Site))
ordihull(benthic.mds, meta_table$Host_spp,col=1:3)
ordispider(benthic.mds, meta_table$Host_spp, col= 1:3, label = TRUE)


#ordiellipse(benthic.mds, Site, col=1:3, kind = "ehull")
#ordispider(benthic.mds, Site, col=1,2,5, label = TRUE)


#### Previous ordination methods


# Dissimilarity martix (transform data)

dist  <- metaMDSiter(benthic.mds, k = 3)
plot(dist)

#

distfa  <- metaMDSiter(fauna_matrix, k = 3) # MDS fauna
plot(distfa)

distfl <- metaMDSiter(flora_matrix, k = 3) # MDS fauna
plot(distfa)

# Clustering and Ordination 

csin <- hclust(all_bray, method="single", main = "") # Single linkage cluster: shortest dissimilarities among cluster members in single linkage (Using previous matrix)
plot(csin, hang = -1, xlab="", sub="") # Plot cluster with 'branches' facing up
rect.hclust(csin, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class

ccom <- hclust(all_bray, method="complete") # Complete linkage cluster: the longest possible dissimilarities in complete linkage
plot(ccom, hang=-1, xlab="", sub="", main = "") # Plot complete linkage plot
rect.hclust(ccom, 3) # Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class

caver <- hclust(all_bray, method="aver") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver, hang=-1, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class


caver <- hclust(all_bray, method="ward.D") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver, hang=-1, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class


caver_all <- hclust(all_bray, method="ward.D2") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver_all, hang=-1, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver_all, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class

# Dendo for flora

caver2 <- hclust(flora_matrix, method="aver") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver2, hang=-1, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver2, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class

caver_flora <- hclust(flora_matrix, method="ward.D2") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver_flora, hang=-1, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver_flora, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class

# Dendo for fauna

caver <- hclust(fauna_matrix, method="aver") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver3, hang=-1, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver3, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class

caver_fauna <- hclust(fauna_matrix, method="ward.D2") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver_fauna, hang=-1, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver_fauna, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class


caver_all <- hclust(all_bray, method="ward.D2") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver_all, hang=-1, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver_all, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class

caver_flora <- hclust(flora_matrix, method="ward.D2") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver_flora, hang=-1, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver_flora, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class

caver_fauna <- hclust(fauna_matrix, method="ward.D2") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver_fauna, hang=-1, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver_fauna, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class

