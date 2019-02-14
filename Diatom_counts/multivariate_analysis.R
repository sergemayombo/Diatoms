# multivariate_analysis.R

library(vegan)
library(BiodiversityR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(tibble)

# https://www.rdocumentation.org/packages/BiodiversityR/versions/2.10-1/topics/nested.anova.dbrda

# Reviewer 1
# The design of the observational study includes 2 treatments - age (young versus old) and host species (Laminaria versus Ecklonia), 4 replicates (4 primary blades from each combination of host algae and age), and 3 subsamples from each blade (pseudoreplicates, if treated incorrectly as replicates). The experimental design is analogous to a 2-way ANOVA, but with community data instead of a single individual response variable. This design can evaluate interactive effects between the two treatments (age and species). The authors’ experimental design is most suited to analyses using PERMANOVA, which is the community statistics version of the ANOVA.

# Please indicate for the readers why the data were transformed and standardised using the stated procedures. Definitely a good idea to transform data, but the readers need to understand why particular procedures were employed. Please describe the Wisconsin double standardisation (row/column standardised by row/column total – to produce relative abundance to total and column/row standardised by column/row max – to produce abundance relative to species max abundance). Why a double standardisation + square-root transformation, as opposed to a single row/column standardization by row/column total + square-root transformation?

# Please indicate for the readers why the data were transformed and standardised using the stated procedures. Definitely a good idea to transform data, but the readers need to understand why particular procedures were employed. Please describe the Wisconsin double standardisation:
# * row/column standardised by row/column total–to produce relative abundance to total and column/row standardised; vs.
# * column/row max–to produce abundance relative to species max abundance.
# Why a double standardisation + square-root transformation, as opposed to a single row/column standardisation by row/column total + square-root transformation?

# AJS: About ANOSIM and PERMANOVA
# "Overall, ANOSIM and the Mantel test were very sensitive to heterogeneity in dispersions, with ANOSIM generally being more sensitive than the Mantel test. In contrast, PERMANOVA and Pillai’s trace were largely unaffected by heterogeneity for balanced designs. [...]. PERMANOVA was also unaffected by differences in correlation structure. [...] PERMANOVA was generally, but not always, more powerful than the others to detect changes in community structure"
# Anderson, M. J., & Walsh, D. C. I. (2013). PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? Ecological Monographs, 83(4), 557–574. http://doi.org/10.1890/12-2010.1

# AJS: About data transformation... Useful when the range of data values is very large. Data are square root transformed, and then submitted to Wisconsin double standardization, or species divided by their maxima, and stands standardized to equal totals. These two standardizations often improve the quality of ordinations, but we forgot to think about them in the initial analysis.


# Get the data in ---------------------------------------------------------

spp <- read.csv(file = "PB_data_matrix.csv", row.names = "Replicate", sep = ",", header = TRUE)
# remove ".spp" from name
colnames(spp) <- str_replace(colnames(spp), "\\.spp", "")
# Logarithmic transformation as suggested by Anderson et al. (2006): log_b (x) + 1 for x > 0, where b is the base of the logarithm; zeros are left as zeros. Higher bases give less weight to quantities and more to presences.
spp.log <- decostand(spp, method = "log")
spp.log.dis <- vegdist(spp.log, method = "bray")
env <- as.tibble(read.csv(file = "PB_diat_env.csv",
                          sep = ",", header = TRUE))
env$plant <- as.factor(env$plant)
env$rep <- as.factor(env$rep)


# betadisper --------------------------------------------------------------

# Before doing the PERMANOVA, first check to see if the dispersion is the same
# Homogeneity of groups
# betadisper studies the differences in group homogeneities
# analogous to Levene’s test of the equality of variances
# can only use one factor as an independent variable
# par(mfrow = c(2, 2))
(mod.spp <- with(env, betadisper(spp.log.dis, host_spp)))
# plot(mod.spp, sub = NULL)
# boxplot(mod.spp)
anova(mod.spp)
permutest(mod.spp)

(mod.size <- with(env, betadisper(spp.log.dis, host_size)))
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
(perm.1 <- adonis2(spp.log.dis ~ (host_spp*host_size)/plant,
                   strata = plant,
                   method = p, data = env))

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
