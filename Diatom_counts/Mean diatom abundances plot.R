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



ggplot(ChickWeight, aes(x = Time, y = weight)) +
  geom_line(aes(group = Chick)) +
  labs(x = expression("cells/mm"^2))

# nMDS

library(vegan)
library(grid)

abund_table<-read.csv("PB_data_matrix.csv",row.names=1,check.names=FALSE)
meta_table<-read.csv("PB_diat_env.csv",row.names=1,check.names=FALSE)

#Just a check to ensure that the samples in meta_table are in the same order as in abund_table
meta_table<-meta_table[rownames(abund_table),]

#Get grouping information
grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
head(grouping_info)

#Get MDS stats
sol<-metaMDS(abund_table,distance = "bray", k = 2, trymax = 50)

#Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
NMDS=data.frame(x=sol$point[,1],y=sol$point[,2],Host_spp=as.factor(grouping_info[,1]),Host_size=as.factor(grouping_info[,2]),Replicate=as.factor(grouping_info[,3]))

#Get spread of points based on Host-spp
plot.new()
ord<-ordiellipse(sol, as.factor(grouping_info[,1]) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()

#Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#Generate ellipse points
df_ell <- data.frame()
for(g in levels(NMDS$Host_spp)){
  if(g!="" && (g %in% names(ord))){
    
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Host_spp==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                  ,Host_spp=g))
  }
}

head(df_ell)
tail(df_ell)

#Generate mean values from NMDS plot grouped on Host_spp
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$Host_spp),mean)

NMDS.mean

#Now do the actual plotting
library(ggplot2)

shape_values<-seq(1,2)

p<-ggplot(data=NMDS,aes(x,y,colour=Host_spp))
p<-p+ annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)
p<-p+ geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)
p<-p+geom_point(aes(shape=Host_size))+scale_shape_manual(values=shape_values)+theme_bw() 
pdf("NMDS.pdf")
print(p)
dev.off()
