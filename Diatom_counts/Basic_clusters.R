# Epiphytic diatoms associated with the south African kelp
# Diatom data exploration, analysis and presentation
# Mayombo Ntamnwe

# 11th Febraury 2018

library(tidyverse)
library(vegan)
library(cluster)
library(ggplot2)
library(ggdendro)
library(tcltk)
library(BiodiversityR)
library(readr)

# read in the data
counts <- read_csv("Diatom_counts_tidy.csv")

# select only some columns
counts.spec <- counts %>%
  select(Site, Host, Replicate, Host_spp, Host_size, Genus, Density) %>%
  na.omit()
counts.spec1 <- Diatom_counts_tidy %>%
  select(Site, Host, Replicate, Host_spp, Host_size, Species, Density) %/%
  na.omit()
counts.gen <- counts %>%
  select(Site, Host, Replicate, Host_spp, Host_size, Genus, Density) %>%
  na.omit()
counts.dens <- counts %>%
  select(Site, Host, Replicate, Host_spp, Host_size, Species, Density) %>%
  na.omit()

summary(counts.spec)

# make into a wide data frame
counts.spec <- spread(counts.spec, key = Genus, value = Density, fill = 0)

# select only some columns
counts <- counts %>%
  select(Site, Host, Replicate, Host_spp, Host_size, Species, Genus, Density) %>%
  na.omit()

diat_counts <- Diatom_counts_tidy %>%
  select(Site, Host, Replicate, Host_spp, Host_size, Species, Genus, Density) %>%
  na.omit()

# make into a wide data frame
counts <- spread(counts, key = Species, value = Density, fill = 0)

# presence/absence only
counts.spec.bin <- decostand(counts.spec[, 7:41], method = "pa")

counts.spec.dis1 <- vegdist(counts.spec.bin, method = "bray", binary = TRUE)

counts.spec.clst1 <- hclust(counts.spec.dis1, method = "ward.D2")

par(mfrow = c(2, 1))
plot(counts.spec.clst1, labels = counts$Host_spp, hang = -1)
plot(counts.clst1, labels = counts$Host_size, hang = -1)
par(mfrow = c(1, 1))

# Bray-Curtis with cell densities
counts.dis2 <- vegdist(counts[, 7:41], method = "bray")

counts.clst2 <- hclust(counts.dis2, method = "ward.D2")

par(mfrow = c(2, 2))
# presence/absence only
plot(counts.clst1, labels = counts$Host_spp, hang = -1, ann = TRUE, xlab = "Host species",
     main = "Presence/absence")
plot(counts.clst1, labels = counts$Host_size, hang = -1, ann = TRUE, xlab = "Host age",
     main = "Presence/absence")
# Bray-Curtis (densities)
plot(counts.clst2, labels = counts$Host_spp, hang = -1, ann = TRUE, xlab = "Host species",
     main = "Cell density")
plot(counts.clst2, labels = counts$Host_size, hang = -1, ann = TRUE, xlab = "Host age",
     main = "Cell density")
par(mfrow = c(1, 2))

# More complex and imaginative analyses are possible, as well as ordination if desired


# Analyses of diatom community strtuctures on the South African kelps

# Shannon and Simpson diversity index based on presence/absence data

head(counts.bin)
tail(counts.bin)
names(counts.bin)
ncol(counts.bin)
nrow(counts.bin)

# Shannon diversity index
shann <- diversity(counts.bin)

# Simpson diveristy index
simp <- diversity(counts.bin, "simpson")

par(mfrow = c(2,2))
hist(shann)
hist(simp)

# Pair-wise distance mesures between samples based on presence/absence data

bray = vegdist(counts.bin, "bray")
gower = vegdist(counts.bin, "gower")
hist(bray)
hist(gower)
par(mfrow = c(1,2))

# Shannon and Simpson diversity index based on abundace data

counts.spec.abund <- counts.spec[, 7:41]
head(counts.spec.abund)
tail(counts.abund)
names(counts.abund)
ncol(counts.abund)

nrow(counts.abund)

# Shannon diversity index
shann_abund <- diversity(counts.abund)

# Simpson diveristy index
simp_abund <- diversity(counts.abund, "simpson")

par(mfrow = c(2,2))
hist(shann_abund)
hist(simp_abund)

# Pair-wise distance mesures between samples based on presence/absence data

bray_abund = vegdist(counts.abund, "bray")

gower_abund = vegdist(counts.abund, "gower")
hist(bray_abund)
hist(gower_abund)
par(mfrow = c(1,2))

# Rarefaction (Rarefy and rarecurve functions) based on species abundance data

specnumber(counts.abund)

sp.abund_1 <- rowSums(counts.abund)
raremax_1 <- min(rowSums(counts.abund))
raremax_1


range(rowSums(counts.abund))
rowSums(counts.abund)
Srare_1 <- rarefy(counts.abund, raremax_1)
par(mfrow = c(1,2))
plot(sp.abund_1, Srare_1, xlab = "Observed No. of species", ylab = "Rarefied No. of species")
abline(0, 1)
rarecurve(counts.abund, step = 20, col = "Blue", cex = 0.6)
          
# Species accumulation curve

rowSums(counts.spec.abund)

par(mfrow = c(1,2))
diat_sp.acc = specaccum(counts.spec.abund, method = "rarefaction")
names(diat_sp.acc)
plot(diat_sp.acc, xvar = "individual", main = "individual based accumulator")

plot(diat_sp.acc, ci.type = "polygon",xlab = "Replicate", main = "confidence polygon", ci.col = "gray50")


diat_sp.acc1 = specaccum(counts.spec.bin, method = "rarefaction")
names(diat_sp.acc1)
par(mfrow = c(1,2))
plot(diat_sp.acc1, xvar = "individual", main = "individual based accumulator")
plot(diat_sp.acc1, ci.type = "polygon", main = "confidence polygon", xlab = "Replicate", ci.col = "gray50")

# Fit non-linear model to species accumulation curves

diat_sp.acc_random = specaccum(counts.spec.bin, method = "random")
diat_sp.acc_nlm = fitspecaccum(diat_sp.acc_random, model = "arrhenius")

names(diat_sp.acc_nlm)
par(mfrow = c(1,1))
plot(diat_sp.acc_nlm, xlab = "Replicate", col = "gray70")
boxplot(diat_sp.acc_random, add = TRUE, xlab = "Replicate", main = "Fit non-linear model to diatom taxa accumulation curves", pch = "+", col = "gray80")

# COmparison between species area curves for subsets of community data

# Species accumulation model on Ecklonia maxima

specaccum(counts.abund[counts$Host == "E_max_A",])

accumresult(counts.abund, y = counts, factor = "Host", level = "E_max_A")

accumresult(counts.abund, y = counts, factor = "Host", level = "E_max_J")

accumresult(counts.abund, y = counts, factor = "Host", level = "L_pal_A")

accumresult(counts.abund, y = counts, factor = "Host", level = "L_pal_J")
accumcomp(counts.spec.abund, y = counts.spec, factor = "Host", method = "exact", conditioned = TRUE)

accumcomp(counts.spec.abund, y = counts.spec, factor = "Host", xlim = c(0, 7), plotit = T)

?accumcomp
dim(counts.abund)
dim(counts)

# Species richness and 95% confidence intervals for kelp associated diatom assemblagesusing four incidence-based estimators

specpool(counts, pool = counts$Host)
(diat_pool_counts = poolaccum(counts.abund))
plot(diat_pool_counts)

# Plotting with ggplot2

library(grid)
library(gridExtra)

ggplot(data = diat_counts, aes(x = diat_counts$Host, y = diat_counts$Density, colour = diat_counts$Diatom_genus))+
  geom_point()

# Comparing subsets of my data for estimated species richness

diat_counts$index = 1:length(diat_counts$Host)
diat_counts.index = as.list(unstack(diat_counts, form = index ~ Host))
diat_counts.index

pacc = function(x, data,...) {poolaccum(data[x,])} 

diat_counts.sp = lapply(diat_counts.index, FUN = pacc, data = diat_counts)
diat_counts.sp

diat_counts.sp$E_max_A
par(mfrow = c(2,2))
plot(diat_counts.sp$E_max_A)
plot(diat_counts.sp$E_max_J)
plot(diat_counts.sp$L_pal_A)
plot(diat_counts.sp$L_pal_J)
par(mfrow = c(1,2))

# Abundance-based richness estimation
eacc = function(x, data,...) {estaccumR(data[x,])}
diat_counts.spe = lapply(diat_counts, FUN = eacc, data = diat_counts)

# Ordination

# Detrended correspondance analysis

library(vegan)
ord <- decorana(counts.spec[, 7:41])
ord
summary(ord)

# Non-metric multidimensional scaling (NMDS)

ord1 <- metaMDS(counts.spec[, 7:41])

ord1

plot(ord1, type = "n")
points(ord1, display = "sites", cex = 1.2, pch = 21, col = "gray50", bg = "gray70")

ggplot() +
  geom_point(data = counts.spec.mds1, aes(x = NMDS1, y = NMDS2, fill = Host), pch = 21, size = 3, colour = NA)

# Analysis of similarity (ANOSIM)
# Host species
counts.ano <- with(counts, anosim(counts.dis1, Host_spp))
plot(counts.ano)
counts.ano
summary(counts.ano)
# Host size
counts.ano1 <- with(counts, anosim(counts.dis1, Host_size))
plot(counts.ano1)
counts.ano1
summary(counts.ano1)
# Host

counts.ano2 <- with(counts, anosim(counts.dis1, Host))
plot(counts.ano2)
counts.ano2
summary(counts.ano2)
# Similarity percentages (SIMPER)
# Host species
sim <- with(counts, simper(counts.abund, Host_spp))
summary(sim)
sim
# Host age
sim1 <- with(counts, simper(counts.abund, Host_size))
summary(sim1)
# Host
sim2 <- with(counts, simper(counts.abund, Host))

summary(sim2)

sim2

# Plotting species abundances
?ggplot

ggplot(data = Diatom_counts_tidy, aes(x = Host, y = Density, fill = Genus)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Host", y = "Density [cells/mm^2]")
  
ggplot(data = diat_counts, aes(x = Host, y = Density, fill = Diatom_genus)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Host", y = "Density (mm-2)")

ggplot(counts.spec, aes(Genera, Density, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_grid(.~Host, drop = TRUE, scales = "free", space = "free_x") +
  theme_bw() + ylab("Density") + xlab("Samples") +
  scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines")) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

ggplot(counts.av, aes(Genus, Density, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_grid(.~Host, drop = TRUE, scales = "free", space = "free_x") +
  theme_bw() + ylab("Diatom density (cells/square millimeter)") + xlab("Diatom genus") +
  scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines")) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

ggplot_alternative <- function()

# Summary stats:
library(Rmisc)
library(ggplot2)
counts.av <- summarySE(counts.spec, measurevar = "Density", groupvars = c("Host", "Genus"), na.rm = TRUE)
counts.av



# Plotting mean diatom abundances with error bars


ggplot(counts.av, aes(Genus, Density, fill = Genus)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Density-se, ymax = Density+se), size = .3, width = .2, position = position_dodge(.9)) +
  facet_grid(.~Host, drop = TRUE, scales = "free", space = "free_x") +
  theme_bw() + ylab("Diatom density (cells/square millimeter)") + xlab("Diatom genus") +
  scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines")) +
  scale_fill_hue(name = "Diatom genus") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
 

ggsave(
  "ggtest.png",
  ggplot_alternative(),
  width = 18,
  height = 15,
  units = "cm",
  dpi = 600
)


# Ordination: Basic method
# Non-metric multidimensional sclaing
library(vegan)
library(MASS)
install.packages()
counts.spec.bin <- decostand(counts.spec[, 7:41], method = "pa")
counts.spec.abund <- counts.spec[, 7:41]

# NMDS with aabundance data
counts.spec.mds1 <- metaMDS(counts.spec.abund, distance = "euclidean", k = 3, autotransform = TRUE)
names(counts.spec.mds1)
counts.spec.mds1

# NMDS with p/a data
counts.spec.mds2 <- metaMDS(counts.spec.bin, distance = "euclidean", k = 3, autotransform = TRUE)
names(counts.spec.mds2)
counts.spec.mds2
# NMDS plot : Sites/samples are shown by black circles, the taxa by red crosses
par(mfrow = c(1,2))
plot(counts.spec.mds1)
plot(counts.spec.mds2)

par(mfrow = c(1,2))
ordiplot(counts.spec.mds1, type = "t")
ordiplot(counts.spec.mds2, type = "t")

plot(counts.spec.mds1, type = "n")
points(counts.spec.mds1, display = "sites", cex = 1.2, pch = 21, col = "gray50", bg = "gray70")

plot(counts.spec.mds2, type = "n")
points(counts.spec.mds2, display = "sites", cex = 1.2, pch = 21, col = "gray50", bg = "gray70")

ggplot(counts.spec.mds1) + 
  geom_point(aes(x = NMDS1, y = NMDS2, col = Replicate, Shape = Host))

library(grid)
counts.keep <- as.numeric(unlist(strsplit(counts[, 8]), ","))

counts.fit <- envfit(counts.spec.mds1, counts[ , 2, drop = F], perm = 999, na.rm = TRUE)
counts.fit

df <- scores(counts.spec.mds1, display = c("sites"))

ggplot() +
  geom_point(data = counts.spec.mds1, aes(NMDS1, NMDS2, colour = "Host"))
plot(df)

?split

# PCA
counts.spec.pca <- rda(counts.spec.abund)
counts.spec.pca

plot(counts.spec.pca)
sum(apply(counts.spec.pca, 2, var))
biplot(counts.spec.pca, scaling = -1)

citation()
