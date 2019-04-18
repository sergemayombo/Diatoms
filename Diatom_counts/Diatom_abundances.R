# Diatom abundances data
# Full dataset
# Data analyses
# Serge Mayombo
# 12 April 2019

# loading libraries
library(readr)
library(lattice)
library(tidyverse)
library(Rmisc)


# Loading datasets
# loading the diatom abundances dataframe

# dat_1 <- read_csv2("Dataset_full.csv")
# dat_1

spp <- read_delim("Dataset_full.csv", 
                           ";", escape_double = FALSE, col_types = cols(`Achnanthes spp` = col_number(), 
                                                                        `Amphora spp` = col_number(), `Asteromphalus spp` = col_number(), 
                                                                        `Cocconeis spp` = col_number(), `Craspedostauros spp` = col_number(), 
                                                                        `Cyclotella spp` = col_number(), 
                                                                        `Cylindrotheca spp` = col_number(), 
                                                                        `Diploneis spp` = col_number(), `Gomphoseptatum spp` = col_number(), 
                                                                        `Grammatophora spp` = col_number(), 
                                                                        `Haslea spp` = col_number(), `Licmophora spp` = col_number(), 
                                                                        `Nagumoea spp` = col_number(), `Navicula spp` = col_number(), 
                                                                        `Nitzschia spp` = col_number(), `Parlibellus spp` = col_number(), 
                                                                        `Rhoicosphenia spp` = col_number(), 
                                                                        `Skeletonema spp` = col_number(), 
                                                                        `Tabularia spp` = col_number(), `Thalassionema spp` = col_number(), 
                                                                        `Trachyneis spp` = col_number()), 
                           trim_ws = TRUE)
View(spp)
# remove ".spp" from name
colnames(spp) <- str_replace(colnames(spp), "\\.spp", "")

# loading the kelp specimen dataframe
env <- read_delim("Dataset_env_full.csv", 
                               ";", escape_double = FALSE, trim_ws = TRUE)

View(env)

# Basic data manupalation and Summary statistics

dim(spp)

dim(env)

summary(spp)

# Data wrangling (Tidy data)
# Merging dataframes
# If we want to retain everything in the Dataset_full_env data frame and merge only what 
# matches in the Dataset_full right data frame, we specify all.x=TRUE. This is known as a
# LEFT JOIN.

diat <- merge(env, spp, all.x=TRUE) # merge the spp and env dataframes 

View(diat)

colnames(spp)

# remove ".spp" from name
colnames(spp) <- str_replace(colnames(spp), "\\.spp", "")

# Converting the merged dataframe into a long dataframe format using function gather_()
keycol <- "diatom_genus" # Name of new key column (made from names of diatom genus)
valuecol <- "Abundance" # Name of new value column = Diatom density
gathercols <- c("Achnanthes spp", "Amphora spp", "Asteromphalus spp", "Cocconeis spp", "Craspedostauros spp", "Cyclotella spp",
                "Cylindrotheca spp", "Diploneis spp", "Gomphoseptatum spp", "Grammatophora spp", "Haslea spp", 
                "Licmophora spp", "Nagumoea spp", "Navicula spp", "Nitzschia spp", "Parlibellus spp", "Rhoicosphenia spp", 
                "Skeletonema spp", "Tabularia spp", "Thalassionema spp", "Trachyneis spp")

diat_long <- gather_(diat, keycol, valuecol, gathercols)

#diat_tidy <- diat %>%
 # gather(Achnanthes spp, Amphora spp, Asteromphalus spp, Cocconeis spp, Craspedostauros spp, Cyclotella spp,
  #        Cylindrotheca spp, Diploneis spp, Gomphoseptatum spp, Grammatophora spp, Haslea spp, Licmophora spp,
   #      Nagumoea spp, Navicula spp, Nitzschia spp, Parlibellus spp, Rhoicosphenia spp, Skeletonema spp, Tabularia spp, 
    #    Thalassionema spp, Trachyneis spp, key = "Genus", value = "Density")

# Summarising the data using the function summarySE
# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval

spp_aver <- summarySE(diat_long, measurevar = "Abundance", groupvars = c("Host_spp", "Host_size", "Part", "diatom_genus"), na.rm = TRUE)

spp_aver

# Plotting mean diatom abundances with error bars representing the standard error of the mean

ggplot(spp_aver, aes(diatom_genus, Abundance, fill = diatom_genus)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Abundance - se, ymax = Abundance + se), size = .3, width = .2, position = position_dodge(.9)) +
  facet_grid(.~Part, drop = TRUE, scales = "free", space = "free_x") +
  theme_bw() + ylab(expression("Diatom density" ~ "(cells mm"^-2*")")) + xlab("Diatom genus") +
  scale_y_continuous(expand = c(0,0)) + theme(strip.background = element_rect(fill="gray85")) +
  theme(panel.margin = unit(0.3, "lines")) +
  scale_fill_hue(name = "Diatom genus", guide = FALSE) +
  theme(axis.text.x = element_text(face = "bold.italic", angle = 90, hjust = 1, vjust = 0.5), strip.text = element_text(size = 6, face = "bold"))
# ggsave("diatom_boxlot1.pdf", width = 6, height = 4, scale = 1.4)

# FIlter observations
# Filter Adult Laminaria observations

ALp <- diat_long %>%
  filter(Host_spp == "Lp", Host_size == "A") # Filter Adult Laminaria pallida observations

ALp

# Summarise Laminaria pallida observations
ALp_aver <- summarySE(ALp, measurevar = "Abundance", groupvars = c("Part", "diatom_genus"), na.rm = TRUE)

ALp_aver

# Plotting mean diatom abundances with error bars representing the standard error of the mean on adult Laminaria pallida

ALp_plot <- ggplot(ALp_aver, aes(diatom_genus, Abundance, fill = diatom_genus)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Abundance - se, ymax = Abundance + se), size = .3, width = .2, position = position_dodge(.9)) +
  facet_grid(.~Part, drop = TRUE, scales = "free", space = "free_x") +
  theme_bw() + ylab(expression("Diatom density" ~ "(cells mm"^-2*")")) + xlab("Diatom genus") +
  scale_y_continuous(expand = c(0,0)) + theme(strip.background = element_rect(fill="gray85")) +
  theme(panel.spacing = unit(0.3, "lines")) +
  scale_fill_hue(name = "Diatom genus", guide = FALSE) +
  theme(axis.text.x = element_text(face = "bold.italic", angle = 90, hjust = 1, vjust = 0.5), strip.text = element_text(size = 6, face = "bold"))

# Filter Juvenile Laminaria pallida
JLp <- diat_long %>%
  filter(Host_spp == "Lp", Host_size == "J") # Filter Juvenile Laminaria pallida observations

JLp

# summarise Jucvenile LP observations 

JLp_aver <- summarySE(JLp, measurevar = "Abundance", groupvars = c("Part", "diatom_genus"), na.rm = TRUE)

JLp_aver

# Plot mean diatom abundance on Juvenile LP with error bars as standard error of the mean
JLp_plot <- ggplot(JLp_aver, aes(diatom_genus, Abundance, fill = diatom_genus)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Abundance - se, ymax = Abundance + se), size = .3, width = .2, position = position_dodge(.9)) +
  facet_grid(.~Part, drop = TRUE, scales = "free", space = "free_x") +
  theme_bw() + ylab(expression("Diatom density" ~ "(cells mm"^-2*")")) + xlab("Diatom genus") +
  scale_y_continuous(expand = c(0,0)) + theme(strip.background = element_rect(fill="gray85")) +
  theme(panel.spacing = unit(0.3, "lines")) +
  scale_fill_hue(name = "Diatom genus", guide = FALSE) +
  theme(axis.text.x = element_text(face = "bold.italic", angle = 90, hjust = 1, vjust = 0.5), strip.text = element_text(size = 6, face = "bold"))

# FIlter Adult Ecklonia maxima

AEm <- diat_long %>%
  filter(Host_spp == "Em", Host_size == "A") # Filter Juvenile Laminaria pallida observations

AEm
# summarise Adult Ecklonia maxima observations 

AEm_aver <- summarySE(AEm, measurevar = "Abundance", groupvars = c("Part", "diatom_genus"), na.rm = TRUE)

AEm_aver

# Plot mean diatom abundance on adult Ecklonia maximum with error bars as standard error of the mean
AEm_plot <- ggplot(AEm_aver, aes(diatom_genus, Abundance, fill = diatom_genus)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Abundance - se, ymax = Abundance + se), size = .3, width = .2, position = position_dodge(.9)) +
  facet_grid(.~Part, drop = TRUE, scales = "free", space = "free_x") +
  theme_bw() + ylab(expression("Diatom density" ~ "(cells mm"^-2*")")) + xlab("Diatom genus") +
  scale_y_continuous(expand = c(0,0)) + theme(strip.background = element_rect(fill="gray85")) +
  theme(panel.spacing = unit(0.3, "lines")) +
  scale_fill_hue(name = "Diatom genus", guide = FALSE) +
  theme(axis.text.x = element_text(face = "bold.italic", angle = 90, hjust = 1, vjust = 0.5), strip.text = element_text(size = 6, face = "bold"))

# FIlter Adult Ecklonia maxima

JEm <- diat_long %>%
  filter(Host_spp == "Em", Host_size == "J") # Filter Juvenile Laminaria pallida observations

JEm
# summarise Adult Ecklonia maxima observations 

JEm_aver <- summarySE(JEm, measurevar = "Abundance", groupvars = c("Part", "diatom_genus"), na.rm = TRUE)

JEm_aver

# Plot mean diatom abundance on adult Ecklonia maximum with error bars as standard error of the mean
JEm_plot <- ggplot(JEm_aver, aes(diatom_genus, Abundance, fill = diatom_genus)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Abundance - se, ymax = Abundance + se), size = .3, width = .2, position = position_dodge(.9)) +
  facet_grid(.~Part, drop = TRUE, scales = "free", space = "free_x") +
  theme_bw() + ylab(expression("Diatom density" ~ "(cells mm"^-2*")")) + xlab("Diatom genus") +
  scale_y_continuous(expand = c(0,0)) + theme(strip.background = element_rect(fill="gray85")) +
  theme(panel.spacing = unit(0.3, "lines")) +
  scale_fill_hue(name = "Diatom genus", guide = FALSE) +
  theme(axis.text.x = element_text(face = "bold.italic", angle = 90, hjust = 1, vjust = 0.5), strip.text = element_text(size = 6, face = "bold"))

# Griding the plots

library(tidyverse)
library(ggpubr)
ggarrange(ALp_plot, JLp_plot, AEm_plot, JEm_plot,
          ncol = 2, nrow = 2, # Set number of rows and columns
          labels = c("A", "B", "C", "D"), # Label each figure
          common.legend = TRUE) # Create common legend

# multivariate_analysis.R

library(permute)
library(vegan)
library(tcltk)
library(BiodiversityR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(tibble)

# standardization method for community ecology with function decostand()
# Logarithmic transformation as suggested by Anderson et al. (2006): log_b (x) + 1 for x > 0, 
# where b is the base of the logarithm; zeros are left as zeros. 
# Higher bases give less weight to quantities and more to presences.

spp.log <- decostand(spp, method = "log", na.rm = FALSE)
spp.log.dis <- vegdist(spp.log, method = "bray")
env_1 <- as.tibble(read.csv(file = "Dataset_env_full.csv",
                          sep = ",", header = TRUE))
# env$plant <- as.factor(env$plant)
env_1$plant <- as.factor(env_1$plant)

# betadisper --------------------------------------------------------------

# Before doing the PERMANOVA, first check to see if the dispersion is the same
# Homogeneity of groups
# betadisper studies the differences in group homogeneities
# analogous to Levene’s test of the equality of variances
# can only use one factor as an independent variable
# par(mfrow = c(2, 2))
(mod.spp <- with(env, betadisper(spp.log.dis, Host_spp)))
# plot(mod.spp, sub = NULL)
# boxplot(mod.spp)
anova(mod.spp)
permutest(mod.spp)

(mod.size <- with(env, betadisper(spp.log.dis, Host_size)))
# plot(mod.size)
# boxplot(mod.size)
anova(mod.size)
permutest(mod.size)


# PERMANOVA ---------------------------------------------------------------

# Permutational multivariate analysis of variance using distance matrices
# (Bray-Curtis similarities by default). ANOSIM uses only ranks of Bray-Curtis,
# so the former preserves more information. PERMANOVA allows for variation
# partitioning and allows for more complex designs (multiple factors, nested
# factors, interactions, covariates, etc.)
# adonis2 studies the differences in the group means
# analogous to multivariate analysis of variance

# note that nestedness should be stated in the model formulation
# as well as in the strata
# "If you have a nested error structure, so that you do not want your data be
# shuffled over classes (strata), you should define strata in your permutation"
# -- Jari Oksannen
(perm.1 <- adonis2(spp.log.dis ~ (Host_spp*Host_size)/rep,
                   strata = rep,
                   method = p, data = env_1)

# within subjects effects
# (perm.2 <- adonis2(spp.log.dis ~ host_spp * host_size + plant,
#                    strata = plant,
#                    method = p, data = env))


# nMDS --------------------------------------------------------------------

spp.nmds <- metaMDS(spp.log, k = 2,trymax = 100,
                    distance = "bray", wascores = TRUE)

scores(spp.nmds, display = "species")
scores(spp.nmds, display = "sites")

(ef <- envfit(spp.nmds, env, permu = 999))
plot(spp.nmds, display = "sites", tck = .02, mgp = c(1.8, 0.5, 0))
plot(ef, p.max = 0.1)

# set things up for the panel of plots (Figure 3)
pdf(file = "../paper_1_v2/Figure__3.pdf", width = 8, height = 7)
col <- c("red3", "blue3")
pch <- c(17, 19)
opar <- par()
plt1 <- layout(rbind(c(1, 1, 2, 2, 3, 3),
                     c(4, 4, 4, 5, 5, 5)),
               heights = c(2, 3),
               respect = TRUE)
layout.show(plt1)
par(mar = c(3,3,1,1))
# plot 1
plot(mod.spp, main = NULL,
     tck = .05, mgp = c(1.8, 0.5, 0), col = col, pch = pch,
     sub = NULL)
# plot 2
plot(mod.size, main = NULL,
     tck = .05, mgp = c(1.8, 0.5, 0), col = col, pch = pch,
     sub = NULL)
# plot 3
stressplot(spp.nmds,
           tck = .05, mgp = c(1.8, 0.5, 0))
# plot 4
par(mar = c(3,3,2,1))
plot(spp.nmds, display = "sites", type = "n",
     main = NULL,
     tck = .05, mgp = c(1.8, 0.5, 0),
     xlim = c(-2, 2), ylim = c(-1, 2))
with(env,
     points(spp.nmds, display = "sites", col = col[host_spp],
            pch = pch[host_spp]))
points(spp.nmds, display = "species", pch = 4, col = "pink")
orditorp(spp.nmds, display = "species", cex = 0.8,
         col = "black", air = 0.01)
with(env,
     ordispider(spp.nmds, groups = host_spp,
                label = TRUE,
                col = col))
with(env, ordiellipse(spp.nmds, groups = host_spp,
                      col = col, label = FALSE))
#plot 5
par(mar = c(3,3,2,1))
plot(spp.nmds, display = "sites", type = "n",
     main = NULL,
     tck = .05, mgp = c(1.8, 0.5, 0),
     xlim = c(-2, 2), ylim = c(-1, 2))
with(env,
     points(spp.nmds, display = "sites", col = col[host_size],
            pch = pch[host_size]))
points(spp.nmds, display = "species", pch = 4, col = "pink")
orditorp(spp.nmds, display = "species", cex = 0.8,
         col = "black", air = 0.01)
with(env,
     ordispider(spp.nmds, groups = host_size,
                label = TRUE,
                col = col))
with(env, ordiellipse(spp.nmds, groups = host_size,
                      col = col, label = FALSE))
dev.off()
par(opar)


# mvabund -----------------------------------------------------------------

library(mvabund)
diat_spp <- mvabund(spp)

# look at the spread of the data using the boxplot function
# not used in paper
par(mar = c(2, 10, 2, 2)) # adjusts the margins
boxplot(spp, horizontal = TRUE, las = 2, main = "Abundance")
dev.off()
par(opar)

# check the mean-variance relationship
# it shows that spp with a high mean also have a high var
meanvar.plot(diat_spp)

# Are there differences in the species composition of the diatom spp. sampled?
# Do some of them specialise on particular spp of kelp, while others are more
# generalised?
# Do some occur more on juveniles, while some are on adults, and which ones
# indiscriminately live across age classes?
# Which species?

# plot(diat_spp ~ env$host_spp, cex.axis = 0.8, cex = 0.8, n.vars = 18)
# plot(diat_spp ~ env$host_size, cex.axis = 0.8, cex = 0.8,
#      type = "bx", n.vars = 18)
# replicated in ggplot2 below...

# scale manually for ggplot2 custom plot
log_fun <- function(x) {
  min_x <- min(x[x != 0], na.rm = TRUE)
  a <- log(x)/min_x
  a[which(!is.finite(a))] <- 0
  return(a)
}

plt1 <- spp %>%
  mutate(host_size = env$host_size) %>%
  gather(key = species, value = abund, -host_size) %>%
  as_tibble() %>%
  group_by(species) %>%
  mutate(log.abund = log_fun(abund)) %>%
  ungroup() %>%
  ggplot(aes(x = fct_reorder(species, abund, .fun = mean), y = log.abund)) +
  geom_boxplot(aes(colour = host_size), size = 0.4, outlier.size = 0,
               fill = "grey90") +
  geom_point(aes(colour = host_size, shape = host_size),
             position = position_dodge2(width = 0.8),
             alpha = 0.6, size = 2.5) +
  scale_colour_manual(name = "Age", values = c("red3", "blue3")) +
  scale_shape_manual(name = "Age", values = c(17, 19)) +
  annotate("text", x = 15, y = 3, size = 4.5,
           label = expression(paste(italic("p"), "=0.017"))) +
  annotate("text", x = 14, y = 3, size = 4.5,
           label = expression(paste(italic("p"), "=0.004"))) +
  scale_y_continuous(name = "Log abundance") +
  coord_flip() + theme_bw() +
  theme(panel.grid.major = element_line(linetype = "dashed",
                                        colour = "turquoise", size = 0.2),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, color = "black",
                                   margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
        axis.text.y = element_text(size = 13, color = "black", face = "italic",
                                   margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
        axis.title.x = element_text(size = 14, vjust = 5.75, color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.ticks = element_line(color = "black", size = 0.5))
ggsave(filename = "../paper_1_v2/Figure__4.pdf", plot = plt1,
       width = 6, height = 10, scale = 0.9)

# try with poisson distribution, but residuals no good
# size_mod1 <- manyglm(diat_spp ~ env$host_spp * env$host_size * env$plant, family = "poisson")
# plot(size_mod1)

# use negative binomial
size_mod2 <- manyglm(diat_spp ~ (env$host_spp*env$host_size)/env$plant,
                     family = "negative binomial")
plot(size_mod2) # better residuals...
anova(size_mod2, test = "wald")
out <- anova(size_mod2, p.uni = "adjusted", test = "wald")
out$table
prop.contrib <- data.frame(spp = colnames(out$uni.test),
                           prop = out$uni.test[3, ],
                           row.names = NULL)
prop.contrib <- prop.contrib %>%
  mutate(perc = round((prop / sum(prop)) * 100, 1)) %>%
  arrange(desc(perc)) %>%
  mutate(cum = cumsum(perc))
