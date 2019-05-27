# Multivariate analysis
# Diatom abundances on South African kelp
# Serge Mayombo
# 24th May 2019

# Impporting data into R

library(vegan)
library(BiodiversityR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(tibble)
library(readr)

?as.tibble
spe <- read.csv2 (file = "Dataset_full.csv", row.names = "Replicate", sep = ";", header = TRUE)

env <- as.tibble(read.csv2 (file = "env_kelp.csv", row.names = "Replicate", sep = ";", header = TRUE))

# Remove samples with 0 diatom abundances

spe <- spe[as.logical(rowSums(spe != 0)), ]

#Just a check to ensure that the samples in meta_table (env) are in the same order as in the species abundance_table (spe)

env <- env[rownames(spe),]

# Get grouping information
#grouping_info <- data.frame(row.names=rownames(spe),t(as.data.frame(strsplit(rownames(spe),"_"))))

#head(grouping_info)

# Logarithmic transformation as suggested by Anderson et al. (2006): 
# log_b (x) + 1 for x > 0, where b is the base of the logarithm; zeros are left as zeros. 
# Higher bases give less weight to quantities and more to presences.

spe.log <- decostand(spe, method = "log") # logarithmic transformation
spe.log.dis <- vegdist(spe.log, method = "bray") # Distance matrix


env$Part <- as.factor(env$Part)
env$rep <- as.factor(env$rep)
env$Host_spp <- as.factor(env$Host_spp)
env$Host_size <- as.factor(env$Host_size)


# NMDS.R
# This script finds a non-parameteric monotonic relationship between the dissimilarities in the samples matrix, 
# and plots the location of each site in a low-dimensional space (similar to principle component analysis).
# loading vegan library

#Get MDS stats
spe.nmds <- metaMDS(spe.log, distance = "bray", k = 2, trymax = 100, wascores = TRUE)

#scores(spe.nmds, display = "species")
#scores(spe.nmds, display = "sites")

#(ef <- envfit(spe.nmds, env, permu = 999))
#plot(spe.nmds, display = "sites", tck = .02, mgp = c(1.8, 0.5, 0))
#plot(ef, p.max = 0.1)

# Make a new data frame, and put Host_spp, Host_size, and Part information there, to be useful for coloring, and shape of points

NMDS = data.frame(NMDS1=spe.nmds$point[,1],NMDS2=spe.nmds$point[,2], Host_spp = as.factor(env$Host_spp), 
                  Host_size = as.factor(env$Host_size), Part = as.factor(env$Part))

# Based on host kelp species
# Get spread of points based on sample specimens 

plot.new()
ord<-ordiellipse(spe.nmds, as.factor(env$Host_spp) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()

#Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. 
#This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 300) 
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

#Generate mean values from NMDS plot grouped on host parts
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$Part),mean)

head(NMDS.mean)
#Now do the actual plotting
library(ggplot2)

shape_values<-seq(1,6)

#p<-ggplot(data=NMDS,aes(NMDS1,NMDS1,colour=Host_spp))
#p<-p+ annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)
#p<-p+ geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)
#p<-p+geom_point(aes(shape=Part))+scale_shape_manual(values=shape_values)+theme_bw() 
#pdf("NMDS_Se.pdf")
#print(p)
#dev.off()

dev.new(title="NMDS plot (Host_spp)")

ggplot(data=NMDS,aes(NMDS1,NMDS2,colour = Host_spp)) +
  annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2) +
  geom_point(aes(shape=Part))+scale_shape_manual(values=shape_values)+theme_bw()


# Based on host kelp size
# Get spread of points based on sample specimens 

plot.new()
ord<-ordiellipse(spe.nmds, as.factor(env$Host_size) ,display = "sites", kind ="sd", conf = 0.95, label = T)
dev.off()

#Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. 
#This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 300) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#Generate ellipse points
df_ell <- data.frame()
for(g in levels(NMDS$Host_size)){
  if(g!="" && (g %in% names(ord))){
    
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Host_size==g,],
                                                     veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                  ,Host_size=g))
  }
}

head(df_ell)
#Generate mean values from NMDS plot grouped on host parts
NMDS.mean=aggregate(NMDS[,1:2],list(group=NMDS$Part),mean)

head(NMDS.mean)
#Now do the actual plotting
library(ggplot2)

shape_values<-seq(1,6)

#p<-ggplot(data=NMDS,aes(NMDS1,NMDS1,colour=Host_spp))
#p<-p+ annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4)
#p<-p+ geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2)
#p<-p+geom_point(aes(shape=Part))+scale_shape_manual(values=shape_values)+theme_bw() 
#pdf("NMDS_Se.pdf")
#print(p)
#dev.off()

dev.new(title="NMDS plot (Host_size)")

ggplot(data=NMDS,aes(NMDS1,NMDS2,colour = Host_size)) +
  annotate("text",x=NMDS.mean$x,y=NMDS.mean$y,label=NMDS.mean$group,size=4) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2), size=1, linetype=2) +
  geom_point(aes(shape=Part))+scale_shape_manual(values=shape_values)+theme_bw()


# nMDS --------------------------------------------------------------------

spp.nmds <- metaMDS(spe.log, k = 2,trymax = 100,
                    distance = "bray", wascores = TRUE)

scores(spp.nmds, display = "species")
scores(spp.nmds, display = "sites")


(ef <- envfit(spp.nmds, env, permu = 999))
plot(spp.nmds, display = "sites", tck = .02, mgp = c(1.8, 0.5, 0))
plot(ef, p.max = 0.1)

col <- c("red3", "blue3", "yellow3", "green3", "cyan3", "purple3")
pch <- c(17, 19)
opar <- par()

dev.new(title="NMDS plot")
plot(spp.nmds, display = "sites", type = "n",
     main = NULL,
     tck = .05, mgp = c(1.8, 0.5, 0),
     xlim = c(-1, 9), ylim = c(-1, 4))
with(env,
     points(spp.nmds, display = "sites", col = col[Part],
            pch = pch[Host_spp]))
points(spp.nmds, display = "species", pch = 4, col = "pink")
orditorp(spp.nmds, display = "species", cex = 0.8,
         col = "black", air = 0.01)
with(env,
     ordispider(spp.nmds, groups = Part,
                label = TRUE,
                col = col))
with(env, ordiellipse(spp.nmds, groups = Host_size,
                      col = col[Host_size], label = FALSE))


# TAXAplot.R
# ============================================================
# Tutorial on drawing a taxa plot using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================

# Loading data

spe <- read.csv2 (file = "Dataset_full.csv", row.names = "Replicate", sep = ";", header = TRUE)

env <- as.tibble(read.csv2 (file = "env_kelp.csv", row.names = "Replicate", sep = ";", header = TRUE))

# Remove samples with 0 diatom abundances

spe <- spe[as.logical(rowSums(spe != 0)), ]


#Just a check to ensure that the samples in meta_table are in the same order as in abund_table
env<-env[rownames(spe),]
head(env)

#Apply proportion normalisation
x<-spe/rowSums(spe)
x<-x[,order(colSums(x),decreasing=TRUE)]
head(x)

#Extract list of top N Taxa
N<-22
taxa_list<-colnames(x)[1:N]
#remove "__Unknown__" and add it to others
#taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
N<-length(taxa_list)

#############################################
# Generate a new table
# Data wrangling (Tidy data)
# Merging dataframes
# If we want to retain everything in the env data frame and merge only what 
# matches in the x data frame at the right, we specify all.x=TRUE. This is known as a
# LEFT JOIN.

# df <- merge(env, p, all.x=TRUE, rownames=T) # merge the p and env dataframes 
#df
#summary(df)
#?ggproto

 ########################################################


#Generate a new table with everything added to Others
#new_x<-data.frame(x[,colnames(x) %in% taxa_list],Others=rowSums(x[,!colnames(x) %in% taxa_list]))


#You can change the Type=grouping_info[,1] should you desire any other grouping of panels
df<-NULL
for (i in 1:dim(x)[2]){
  tmp<-data.frame(row.names=NULL,Sample=rownames(x),Taxa=rep(colnames(x)[i],dim(x)[1]),Value=x[,i],Type1=env$Host_spp, Type2=env$Host_size)
  if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5",
             "#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2",
             "#00998F","#740AFF","#990000","#FFFF00");

dev.new(title="Diatom taxa abundances")

ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~Type, drop=TRUE,scale="free",space="free_x") +
  scale_fill_manual(values=colours[1:(N+1)]) +
  theme_bw()+ylab("Proportions") +
  scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+
  theme(panel.margin = unit(0.3, "lines")) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

# Since the dataset is a bit large
# we'll plot diatom abundance proportions on adults and juveniles specimens of both host kelp species separately
# filter()` allows you to subset observations based on their values. The first argument is the name of the data frame. 
# The second and subsequent arguments are the expressions that filter the data frame.
# load library

library(lubridate)

df_ad <- filter(df, Type2 == "A")
df_juv <- filter(df, Type2 == "J")

# Now we can plot diatom abundance proportions on adult and juvenile specimens separately.
# Diatom abundance on adult kelp specimens

dev.new(title="Diatom abundances on adult specimens")

ggplot(df_ad,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type1, drop=TRUE,scale="free",space="free_x")+
  scale_fill_manual(values=colours[1:(N+1)]) + theme_bw()+ylab("Proportions") +
  scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines")) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size = 7), strip.text = element_text(size = 12, face = "bold"))+
  theme(legend.key.size = unit(0.5, "cm"), legend.key.width = unit(0.5,"cm"), legend.text = element_text(size = 8, face = "italic"))

# pdf("Taxaplot1.pdf",height=6,width=21)
# print(p)
# dev.off()

# Now we'll plot diatom abundance proportions on juvenile specimens
dev.new(title="Diatom abundances on juvenile specimens")

ggplot(df_juv,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Type1, drop=TRUE,scale="free",space="free_x")+
  scale_fill_manual(values=colours[1:(N+1)]) + theme_bw()+ylab("Proportions") +
  scale_y_continuous(expand = c(0,0))+theme(strip.background = element_rect(fill="gray85"))+theme(panel.margin = unit(0.3, "lines")) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size = 7), strip.text = element_text(size = 12, face = "bold"))+
  theme(legend.key.size = unit(0.5, "cm"), legend.key.width = unit(0.5,"cm"), legend.text = element_text(size = 8, face = "italic"))

# KW.R ********************************************
# ============================================================
# Tutorial on finding significant taxa and then plotting them using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================
library(reshape)
library(ggplot2)

#Load the abundance table 
#abund_table<-read.csv2("Dataset_full.csv",row.names=1,check.names=FALSE)
abund_table<- read.csv2 (file = "Dataset_full.csv", row.names = "Replicate", sep = ";", header = TRUE)

env <- as.tibble(read.csv2 (file = "env_kelp.csv", row.names = "Replicate", sep = ";", header = TRUE))

#Transpose the data to have sample names on rows
#abund_table<-t(abund_table)

# Remove samples with 0 diatom abundances

abund_table <- abund_table[as.logical(rowSums(abund_table != 0)), ]


#Get grouping information
#grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
# > head(grouping_info)
# X1 X2 X3
# T_2_1   T  2  1
# T_2_10  T  2 10
# T_2_12  T  2 12
# T_2_2   T  2  2
# T_2_3   T  2  3
# T_2_6   T  2  6

#Use countries as grouping information
groups<-as.factor(env$Host_size)

#Apply normalisation (either use relative or log-relative transformation)
#data<-abund_table/rowSums(abund_table)
data<-log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
data<-as.data.frame(data)

#Reference: http://www.bigre.ulb.ac.be/courses/statistics_bioinformatics/practicals/microarrays_berry_2010/berry_feature_selection.html
kruskal.wallis.alpha=0.01
kruskal.wallis.table <- data.frame()
for (i in 1:dim(data)[2]) {
  ks.test <- kruskal.test(data[,i], g=groups)
  # Store the result in the data frame
  kruskal.wallis.table <- rbind(kruskal.wallis.table,
                                data.frame(id=names(data)[i],
                                           p.value=ks.test$p.value
                                ))
  # Report number of values tested
  cat(paste("Kruskal-Wallis test for ",names(data)[i]," ", i, "/", 
            dim(data)[2], "; p-value=", ks.test$p.value,"\n", sep=""))
}


kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]

kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, 
                                    size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)

kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value,
                                                   decreasing=FALSE), ]
kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
pdf("KW_correction3.pdf")
plot(kruskal.wallis.table$p.value,
     kruskal.wallis.table$E.value,
     main='Multitesting corrections',
     xlab='Nominal p-value',
     ylab='Multitesting-corrected statistics',
     log='xy',
     col='blue',
     panel.first=grid(col='#BBBBBB',lty='solid'))
lines(kruskal.wallis.table$p.value,
      kruskal.wallis.table$FWER,
      pch=20,col='darkgreen', type='p'
)
lines(kruskal.wallis.table$p.value,
      kruskal.wallis.table$q.value,
      pch='+',col='darkred', type='p'
)
abline(h=kruskal.wallis.alpha, col='red', lwd=2)
legend('topleft', legend=c('E-value', 'p-value', 'q-value'), col=c('blue', 'darkgreen','darkred'), lwd=2,bg='white',bty='o')
dev.off()

last.significant.element <- max(which(kruskal.wallis.table$q.value <= kruskal.wallis.alpha))
selected <- 1:last.significant.element
diff.cat.factor <- kruskal.wallis.table$id[selected]
diff.cat <- as.vector(diff.cat.factor)

print(kruskal.wallis.table[selected,])

#Now we plot taxa significantly different between the categories
df<-NULL
for(i in diff.cat){
  tmp<-data.frame(data[,i],groups,rep(paste(i," q = ",round(kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"],5),sep=""),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Type","Taxa")

p<-ggplot(df,aes(Type,Value,colour=Type))+ylab("Log-relative normalised")
p<-p+geom_boxplot()+geom_jitter()+theme_bw()+
  facet_wrap( ~ Taxa , scales="free", ncol=3)
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
pdf("KW_significant2.pdf",width=10,height=14)
print(p)
dev.off()

# Cluster.R
# ============================================================
# Tutorial on hierarchical clustering and plotting samples using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================

library(ggplot2)
library(ggdendro)

#Load the abundance table 
#abund_table<-read.csv2("Dataset_full.csv",row.names=1,check.names=FALSE)
abund_table<- read.csv2 (file = "Dataset_full.csv", row.names = "Replicate", sep = ";", header = TRUE)

env <- as.tibble(read.csv2 (file = "env_kelp.csv", row.names = "Replicate", sep = ";", header = TRUE))

#Transpose the data to have sample names on rows
#abund_table<-t(abund_table)

# Remove samples with 0 diatom abundances

abund_table <- abund_table[as.logical(rowSums(abund_table != 0)), ]


#Get grouping information
#grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
# > head(grouping_info)
# X1 X2 X3
# T_2_1   T  2  1
# T_2_10  T  2 10
# T_2_12  T  2 12
# T_2_2   T  2  2
# T_2_3   T  2  3
# T_2_6   T  2  6

#Use countries as grouping information
env$spp<-as.factor(env$Host_spp)
env$size<-as.factor(env$Host_size)
env$Part<-as.factor(env$Part)


#Load the abundance table 
#abund_table<-read.csv("SPE_pitlatrine.csv",row.names=1,check.names=FALSE)

#Transpose the data to have sample names on rows
#abund_table<-t(abund_table)

#Get grouping information
#grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
# > head(grouping_info)
# X1 X2 X3
# T_2_1   T  2  1
# T_2_10  T  2 10
# T_2_12  T  2 12
# T_2_2   T  2  2
# T_2_3   T  2  3
# T_2_6   T  2  6

betad<-vegdist(abund_table,method="bray")

# Use Adonis to test for overall differences
#res_adonis <- adonis(betad ~ env$spp, env$spp) 

#Cluster the samples
hc <- hclust(betad)

#We will color the labels according to kelp parts
hc_d <- dendro_data(as.dendrogram(hc))
hc_d$labels$Type<-env$Part         #[as.character(hc_d$labels$label),1]

#Coloring function
gg_color_hue<-function(n){
  hues=seq(15,375,length=n+1)
  hcl(h=hues,l=65,c=100)[1:n]
}

cols=gg_color_hue(length(unique(hc_d$labels$Type)))
hc_d$labels$color=cols[hc_d$labels$Type]

## Plot clusters
p1 <- ggplot(data = segment(hc_d)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  scale_x_discrete(labels=label(hc_d)$label) +
  ylab("Distance (beta diversity = bray)") + theme_bw()+
  theme(axis.text.y = element_text(color = hc_d$labels$color),
        axis.title.y = element_blank())
p1 <- p1 + geom_point(data=hc_d$label, aes(x = x, y = y, color = Type), inherit.aes =F, alpha = 0)
p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
  scale_color_manual(values = cols)
pdf("Cluster.pdf",height=10)
print(p1)
dev.off()


## Plot clusters

dev.new(title="cluster analysis")
ggplot(data = segment(hc_d)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  coord_flip() +
  scale_x_discrete(labels=label(hc_d)$label) +
  ylab("Distance (beta diversity = bray)") + theme_bw()+
  theme(axis.text.y = element_text(color = hc_d$labels$color),
        axis.title.y = element_blank()) +
  geom_point(data=hc_d$label, aes(x = x, y = y, color = Type), inherit.aes =F, alpha = 0) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
  scale_color_manual(values = cols)

#pdf("Cluster.pdf",height=10)
print(p1)
dev.off()

## Similarity profile analysis (SIMPROF)
library(clustsig)

simprof.plot(hc_d, leafcolors=NA, plot=TRUE, fill=TRUE,
             leaflab="perpendicular", siglinetype=1)

library(vegan)

# calculate the dissimilarity matrix with vegdist so you can use the sorenson/bray 
#method
distBray <- vegdist(abund_table, method = "bray") 

# calculate the clusters with ward.D2
clust1 <- hclust(distBray, method = "ward.D2")
clust1

#plot the cluster dendrogram with dendextend
library(dendextend)
library(ggdendro)
library(ggplot2)

dev.new(title="cluster analysis")
dend <- clust1 %>% as.dendrogram %>%
  set("branches_k_color", k = 5) %>% set("branches_lwd", 0.5)  %>%  set("clear_leaves") %>% 
  set("labels_colors", k = 5)  %>% set("leaves_cex", 0.5) %>%
  set("labels_cex", 0.5)
ggd1 <- as.ggdend(dend)
ggplot(ggd1, horiz = TRUE)

# Run simprof on the data
res <- simprof(data=spe.log, 
               method.distance="actual_braycurtis")
# Graph the result
pl.color <- simprof.plot(res)


## Not run:
# Load the USArrests dataset included with R
# And use abbreviations of state names
# We leave out the third column because
# it is on a different scale
usarrests<-USArrests[,c(1,2,4)]
rownames(usarrests)<-state.abb
# Run simprof on the data
res <- simprof(data=usarrests,
               method.distance="euclidean")
# Graph the result
pl.color <- simprof.plot(res)


# NB.R
# The DESeq {DESeq2} package allows negative binomial GLM fitting and Wald statistics for abundance data. You can use this script as an alternative to KW.R to find taxa that are significantly different between different conditions.
# ============================================================
# Tutorial on finding significant taxa and then plotting them using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)
# =============================================================

library(ggplot2)
library(DESeq2)

#Load the abundance table 
abund_table<- read.csv2 (file = "Dataset_full.csv", row.names = "Replicate", sep = ";", header = TRUE)

env <- as.tibble(read.csv2 (file = "env_kelp.csv", row.names = "Replicate", sep = ";", header = TRUE))

#Transpose the data to have sample names on rows
#abund_table<-t(abund_table)

# Remove samples with 0 diatom abundances

abund_table <- abund_table[as.logical(rowSums(abund_table != 0)), ]
#Get grouping information
grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
# > head(grouping_info)
# X1 X2 X3
# T_2_1   T  2  1
# T_2_10  T  2 10
# T_2_12  T  2 12
# T_2_2   T  2  2
# T_2_3   T  2  3
# T_2_6   T  2  6

#We will convert our table to DESeqDataSet object
countData = round(as(abund_table, "matrix"), digits = 0)
# We will add 1 to the countData otherwise DESeq will fail with the error:
# estimating size factors
# Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
# every gene contains at least one zero, cannot compute log geometric means
countData<-((countData+1)) 

dds <- DESeqDataSetFromMatrix(countData, grouping_info, as.formula(~ X1))

#Reference:https://github.com/MadsAlbertsen/ampvis/blob/master/R/amp_test_species.R

#Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
#Some reason this doesn't work: data_deseq_test = DESeq(dds, test="wald", fitType="parametric")
data_deseq_test = DESeq(dds)

## Extract the results
res = results(data_deseq_test, cooksCutoff = FALSE)
res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]), OTU = rownames(res))

sig = 0.001
fold = 0
plot.point.size = 2
label=T
tax.display = NULL
tax.aggregate = "OTU"

res_tax_sig = subset(res_tax, padj < sig & fold < abs(log2FoldChange))

res_tax_sig <- res_tax_sig[order(res_tax_sig$padj),]

## Plot the data
### MA plot
res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
res_tax$Significant[is.na(res_tax$Significant)] <- "No"
p1 <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) +
  geom_point(size = plot.point.size) +
  scale_x_log10() +
  scale_color_manual(values=c("black", "red")) +
  labs(x = "Mean abundance", y = "Log2 fold change")+theme_bw()
if(label == T){
  if (!is.null(tax.display)){
    rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
  } else {
    rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
  }
  p1 <- p1 + geom_text(data = subset(rlab, Significant == "Yes"), aes(label = Display), size = 4, vjust = 1)
}
pdf("NB_MA.pdf")
print(p1)
dev.off()

res_tax_sig_abund = cbind(as.data.frame(countData[rownames(res_tax_sig), ]), OTU = rownames(res_tax_sig), padj = res_tax[rownames(res_tax_sig),"padj"]) 

#Apply normalisation (either use relative or log-relative transformation)
#data<-abund_table/rowSums(abund_table)
data<-log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
data<-as.data.frame(data)

#Now we plot taxa significantly different between the categories
df<-NULL
for(i in res_tax[rownames(res_tax_sig),"OTU"]){
  tmp<-data.frame(data[,i],groups,rep(paste(i," padj = ",round(res_tax[i,"padj"],5),sep=""),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Type","Taxa")

p<-ggplot(df,aes(Type,Value,colour=Type))+ylab("Log-relative normalised")
p<-p+geom_boxplot()+geom_jitter()+theme_bw()+
  facet_wrap( ~ Taxa , scales="free", ncol=3)
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
pdf("NB_significant.pdf",width=10,height=14)
print(p)
dev.off()

