# Epiphytic diatoms associated with the south African kelp
# Diatom data exploration, analysis and presentation
# Serge Mayombo

# 10th May 2018
library(tidyverse)
library(readr)
## get packages installed
## packs = as.data.frame(installed.packages(.libPaths()[1]), stringsAsFactors = F)

## and now re-install install packages using install.packages()
## install.packages(packs$Package)

# read in the data
counts <- read_csv("Diatom_counts_tidy.csv")

# select only some columns
counts.spec <- counts %>%
  select(Site, Host, Replicate, Host_spp, Host_size, Genus, Density) %>%
  na.omit()

# Summary stats:
library(Rmisc)
library(ggplot2)

counts.av <- summarySE(counts.spec, measurevar = "Density", groupvars = c("Host", "Genus"), na.rm = TRUE)

counts.av



# Plotting mean diatom abundances with error bars

ggplot(counts.av, aes(Genus, Density, fill = Genus)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Density - se, ymax = Density + se), size = .3, width = .2, position = position_dodge(.9)) +
  facet_grid(.~Host, drop = TRUE, scales = "free", space = "free_x") +
  theme_bw() + ylab(expression("Diatom density" ~ "(cells mm"^-2*")")) + xlab("Diatom genus") +
  scale_y_continuous(expand = c(0,0)) + theme(strip.background = element_rect(fill="gray85")) +
  theme(panel.margin = unit(0.3, "lines")) +
  scale_fill_hue(name = "Diatom genus", guide = FALSE) +
  theme(axis.text.x = element_text(face = "bold.italic", angle = 90, hjust = 1, vjust = 0.5), strip.text = element_text(size = 6, face = "bold"))
ggsave("diatom_boxlot1.pdf", width = 6, height = 4, scale = 1.4)


# Non-metric Multidimensional Scaling (nMDS)
#library(ecodist)
library(vegan)

## mds
## Data preparation
abund_data <- read.csv(file = "PB_data_matrix.csv", row.names = 'Replicate', sep = ",", header = TRUE)

#Site  <- abund_data$Replicate
meta_table <- read.csv(file = "PB_diat_env.csv",row.names= "Replicate",sep = ",", header = TRUE)

# Constrained coordination analysis
ord1 <- cca(abund_data ~ Host_spp + Host_size, data = meta_table)
ord1
# ... 20.7% of the variation in the diatom abundance data is explained by the
# constraining variables, Host_spp + Host_size
# ... there are two canonical axes, CCA1 and CCA2 (i.e. the constrained axes)
coef(ord1)


# Plot the centroids for factor constraints and site scores: coordinates of the
# sites as expressed in the space of the response variables Y (i.e. abund_data)...
# Scaling 1 – distance biplot:
plot(ord1, display = c("sp", "cn", "lc"), scaling = 1)

plot(ord1)

# Partitioning the variance--do a permutation tests for the significance of
# constraints; all constraints tested simultaneously:
anova.cca(ord1)

# Do a permutation test for significance of individual terms:
anova.cca(ord1, by = "term", permutations = 199)

# Similar to above:
anova.cca(ord1, by = "mar", permutations = 199)

# Test for significance of axes:
anova(ord1, by = "axis", permutations = 499)
# ...only the first canonical axis is significant; see the figure, plot(ord1)...


# Distance matrix
all_bray <- vegdist(abund_data, method = "bray") # Create bray - curtis dissimilarity matrix

# NMDS
benthic.mds <- metaMDS(abund_data)
benthic.mds

#####################################

data.scores <- as.data.frame(scores(benthic.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(benthic.mds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data


## Plotting mds: Method 1
## Host_spp
pdf("NMDS1.pdf")
plot(benthic.mds, type = "n")
points(benthic.mds, display = "sites", cex = 0.8, pch = 21, col = c(fill = meta_table$Host_spp))
ordihull(benthic.mds, meta_table$Host_spp, col=1:3)
ordiellipse(benthic.mds, meta_table$Host_spp, col=1:3, draw="polygon")
ordispider(benthic.mds, meta_table$Host_spp, col=1:3, label = TRUE)
dev.off()

## Host_size
pdf("NMDS2.pdf")
plot(benthic.mds, type = "n")
points(benthic.mds, display = "sites", cex = 0.8, pch = 21, col = c(fill = meta_table$Host_size))
ordiellipse(benthic.mds, meta_table$Host_size, col=1:3, draw="polygon")
ordispider(benthic.mds, meta_table$Host_size, col=1:3, label = TRUE)
ordihull(benthic.mds, meta_table$Host_size, col=1:3)
dev.off()

# Basic cluster analysis

caver_all <- hclust(all_bray, method="ward.D2") # Average linkage cluster: distances among cluster centroids in average linkage
pdf("hcluster.pdf", width = 6, height = 4)
plot(caver_all, hang=-1, cex = 0.6, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver_all, 3)#
grp <- cutree(caver_all, 3)
dev.off()


## ordination diagram of clustering groups
pdf("ord.pdf")
ord <- cca(abund_data)
plot(ord, display = "sites")
ordihull(ord, grp, lty = 2, col = "red")
dev.off()

# Complete method clustering
ccom <- hclust(all_bray, method="complete") # Complete linkage cluster: the longest possible dissimilarities in complete linkage
pdf("hcclustercom.pdf", width = 6, height = 4) # Save as pdf document
plot(ccom, hang=-1, cex = 0.6, xlab="", sub="", main = "") # Plot complete linkage plot
rect.hclust(ccom, 3) # Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class
dev.off()

## Analysis of similarity (ANOSIM) and Similarity Percentage (SIMPER)
## Data preparation

# $Host_spp
abund.dist <- vegdist(abund_data, method = "bray")
abund.ano <- anosim(abund.dist, meta_table$Host_spp, permutations = 999)
abund.ano
summary(abund.ano)
plot(abund.ano)


# $Host_size
abund.ano1 <- anosim(abund.dist, meta_table$Host_size, permutations = 999)
summary(abund.ano1)
plot(abund.ano1)

## Permuational multivariate analysis of variance
#Host_spp:Host_size
(adonis(abund.dist ~ Host_spp * Host_size, method = "bray", data = meta_table, permutations = 999))

# Simper
# $Host_spp
(sim <- with(meta_table, simper(abund_data, Host_spp), ordered = TRUE, permutations = 999))
summary(sim)

# $Host_size
(sim1 <- with(meta_table, simper(abund_data, Host_size), ordered = T, permutations = 999))
summary(sim1)


#### Minimum, maximum and average temperatures

read.csv("winter_temp.csv")

#######################
# Cluster analysis

taxa.norm <- decostand(abund_data, "normalize")
taxa.ch <- vegdist(taxa.norm, "bray")
taxa.ward <- hclust(taxa.ch, method = "ward.D")
plot(taxa.ward)
taxa.com <- hclust(taxa.ch, method = "complete")
plot(taxa.com)
# grouping
k <- 3
taxa.ward.g <- cutree (taxa.ward, k)
taxa.com.g <- cutree(taxa.com, k)
## Complete linkage vs ward
table(taxa.com.g, taxa.ward.g)
# Reorder dendrogram so that objects' order in the dissimilarity matrix is respected as much as possible
taxa.chwo <- reorder(taxa.ward, taxa.ch)

# plot reordered dendrogram with group labels
plot(taxa.chwo, hang = -1, xlab = "4 groups", sub = "", ylab = "Height", main = "", labels = cutree(taxa.chwo, k=k))
rect.hclust(taxa.chwo, k=k)

# Heat map of the distance matrix ordered with the dendrogram
pdf("heatmap1.pdf")
dend <- as.dendrogram(taxa.chwo)
heatmap(as.matrix(taxa.ch), Rowv = dend, symm = TRUE, margin =c(4,4))
dev.off()

# ordered community table
or <- vegemite(abund_data, taxa.chwo, scale = "log")

# Heat map of diatom abundance community table with dendrogram
#install.packages("RColorBrewer")
library(RColorBrewer)

pdf("heatmapTaxa.pdf")
heatmap(t(abund_data[rev(or$species)]), Rowv = NA, Colv = dend, col = c("white", brewer.pal(9, "Reds")), scale = "none", margin = c(5,10), ylab = "Diatom taxa (weighted averages on samples)", xlab = "Samples")
dev.off()

pdf("diat_abund.pdf")
tabasco(abund_data, taxa.com, xlab = "Samples")
dev.off()

#################################################################

