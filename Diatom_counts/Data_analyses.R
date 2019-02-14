# Epiphytic diatoms associated with the south African kelp
# Diatom data exploration, analysis and presentation
# Serge Mayombo

# 11th Febraury 2018
library(tidyverse)
library(readr)

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
  theme_bw() + ylab("Diatom density (cells/square millimeter)") + xlab("Diatom genus") +
  scale_y_continuous(expand = c(0,0)) + theme(strip.background = element_rect(fill="gray85")) +
  theme(panel.margin = unit(0.3, "lines")) +
  scale_fill_hue(name = "Diatom genus", guide = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.text = element_text(size = 6, face = "bold"))
ggsave("diatom_boxlot.png", width = 6, height = 4, scale = 1.4)


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
library(ecodist)
library(ggplot2)
library(vegan)
library(reshape2)
library(zoo)
library(readr)

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

# Basic cluster analysis

caver_all <- hclust(all_bray, method="ward.D2") # Average linkage cluster: distances among cluster centroids in average linkage
pdf("hcluster.pdf", width = 6, height = 4)
plot(caver_all, hang=-1, xlab="", sub="", main = "") # Plot average linkage cluster
rect.hclust(caver_all, 3)# 
dev.off()
