## Analysis of similarity (ANOSIM) and Similariaty Percentage (SIMPER)
## Data preparation

library(readr)
library(vegan)

abund_data <-  read.csv(file = "PB_data_matrix.csv", row.names = 'Replicate', sep = ",", header = TRUE) 

#Site  <- abund_data$Replicate
meta_table<-read.csv(file = "PB_diat_env.csv",row.names= "Replicate",sep = ",", header = TRUE)
# alldata2 <- subset(alldata2,select=-c(Site))

all_bray <- vegdist(abund_data, method = "bray") # Create bray - curtis dissimilarity matrix
all_bray

# Anosim 
# $Host_spp
abund.dist <- vegdist(abund_data, method="bray")
attach(meta_table)
abund.ano <- anosim(abund.dist, Host_spp, permutations = 999)
abund.ano
summary(abund.ano)
plot(abund.ano)


# $Host_size

abund.ano1 <- anosim(abund.dist, Host_size, permutations = 999)
summary(abund.ano1)
plot(abund.ano1)

#Host_spp:Host_size
(adonis(abund.dist ~ Host_spp*Host_size, method="bray", data=meta_table, permutations=999))
(adonis2(abund.dist ~ Host_spp*Host_size, method="bray", data=meta_table, permutations=999))

# Permutational Multivariate analysis
#Host_size:Host_spp
(adonis(abund.dist ~ Host_size*Host_spp, method="bray", data=meta_table, permutations=999))
adonis1 <- adonis2(abund.dist~Host_size/Host_spp, method="bray", strata = Host_size, data = meta_table)
adonis1

(adonis2(abund_data ~ Host_spp*Host_size, method="bray", data=meta_table, permutations=999))


# Simper
# $Host_spp
(sim <- with(meta_table, simper(abund_data, Host_spp), ordered = TRUE, permutations = 999))
summary(sim)

# $Host_size
(sim1 <- with(meta_table, simper(abund_data, Host_size), ordered=T, permutations = 999))
summary(sim1)
