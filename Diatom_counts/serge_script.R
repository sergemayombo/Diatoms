# ANOSIM 

ano <- anosim(Alldata_wide.1, grouping = Site, permutations = 999, distance = "bray", strata = NULL, parallel = NULL) # ANOSIM for data between sites
summary(ano)
plot(ano) # Plotting results of ANOSIM

# SIMPER
sim <- simper(Alldata_wide.1, group = Site, permutations = 999, trace = FALSE, parallel = 1) # Similarity percentages (SIMPER) for all data across sites NOTE: Make sure parallel setting is correct
sim
summary(sim)

# Transform abundances to presence-absence (1-0)
all_pa <- decostand(alldata_wide, method="pa") 

# cluster analysis

caver_all <- hclust(all_bray, method="ward.D2") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver_all, hang=-1, xlab="", sub="", main = "", cex = 1.2, cex.lab = 1.2) # Plot average linkage cluster
rect.hclust(caver_all, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class
caver_flora <- hclust(flora_matrix, method="ward.D2") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver_flora, hang=-1, xlab="", sub="", main = "", cex = 1.2, cex.lab = 1.2) # Plot average linkage cluster
rect.hclust(caver_flora, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class
caver_fauna <- hclust(fauna_matrix, method="ward.D2") # Average linkage cluster: distances among cluster centroids in average linkage
plot(caver_fauna, hang=-1, xlab="", sub="", main = "", cex = 1.2, cex.lab = 1.2) # Plot average linkage cluster
rect.hclust(caver_fauna, 3)# Hierarchic clustering method: The extremes are all observations in a single class, and each observation in its private class

## ANOVA

shaprio_length <- shapiro.test(kelp_all_chp3$Length)
shaprio_length
capture.output(shaprio_length, file="shaprio_length.doc")
anova_length <- aov(Length ~ Site * Plot, data = kelp_length_chp3)
sum_anova_length <- summary(anova_length)
sum_anova_length
capture.output(sum_anova_length, file="anova_length.doc")
print(anova_length)
TukeyHSD_length <- TukeyHSD(anova_length)
TukeyHSD_length
capture.output(TukeyHSD_length, file = "TukeyHSD_length.doc")
model.tables(anova_length, type = "means", se = TRUE)

kruskal.test(Length ~ Site, data = kelp_length_chp3)
posthoc.kruskal.nemenyi.test(Length ~ Site, data = kelp_length_chp3, method="Tukey")

## Boxplot of adults and points

means_adults <- aggregate(Adults ~ Site, kelp_all, mean) # This calculates the mean so you can add it to graph

box_adults_colour <- ggplot(data = kelp_all, aes(x=Site, y=Adults)) + geom_boxplot(notch = TRUE, outlier.shape=NA,aes(fill=Site), width = 0.3) + theme(axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0.3)) + theme(panel.background = element_rect(fill = 'white', colour = 'grey'))+ xlab("Sites") + ylab("Count of adult kelp") + ggtitle("") + geom_jitter(position = position_jitter(width = 0.3, height = 0.3), alpha = 5/10) + theme(aspect.ratio = 0.8) 
box_adults_colour
ggsave("box_adults_colour.png", width = 18, height = 15, units = "cm", dpi = 100)



