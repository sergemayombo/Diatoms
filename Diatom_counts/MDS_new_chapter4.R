
setwd(dir = "C:/Users/Ross/Documents/R/Work directory/recovery/")

library(vegan)
library(MASS)

## Bordtjies
## Non-metric Multidimensional scaling
## All data
### !!! Chooose a data set first !!!

recovery_all_B <-  read.csv(file = "recovery_all_B_std.csv", row.names = 'ID', sep = ",", header = TRUE) 
Phase  <- recovery_all_B$Phase
recovery_all_B <- subset(recovery_all_B,select=-c(Phase))

recovery_all_O <-  read.csv(file = "recovery_all_O_std.csv", row.names = 'ID', sep = ",", header = TRUE) 
Phase  <- recovery_all_O$Phase
recovery_all_O <- subset(recovery_all_O,select=-c(Phase))


benthic.mds <- metaMDS(recovery_all_O)
benthic.mds

data.scores <- as.data.frame(scores(benthic.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(benthic.mds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

## Plotting mds: chapter4

plot(benthic.mds, type = "p")
ordispider(benthic.mds, Phase, col= 1:4, label = TRUE)
ordihull(benthic.mds, Phase,col=1:4)


multiplot() # Remember to activate function before running


