library(vegan)
library(ggplot2)
#data(dune)
#View(dune)

abund_table<-read.csv("PB_data_matrix.csv",row.names=1,check.names=FALSE)
meta_table<-read.csv("PB_diat_env.csv",row.names=1,check.names=FALSE)

# calculate distance for NMDS
sol <- metaMDS(abund_table, distance = "bray", k = 2, trymax = 50)
sol# Create meta data for grouping
#MyMeta = data.frame(
 # sites = c(2,13,4,16,6,1,8,5,17,15,10,11,9,18,3,20,14,19,12,7),
 # amt = c("hi", "hi", "hi", "md", "lo", "hi", "hi", "lo", "md", "md", "lo", 
  #        "lo", "hi", "lo", "hi", "md", "md", "lo", "hi", "lo"),
  #row.names = "sites")


# plot NMDS using basic plot function and color points by "amt" from MyMeta
plot(sol$points, col = meta_table$Host_spp)

# draw dispersion ellipses around data points
ordiellipse(sol, meta_table$Host_spp, display = "sites", kind = "sd", label = T, col = 1:3)
ordihull(benthic.mds, meta_table$Host_spp,col=1:3)
ordispider(benthic.mds, meta_table$Host_spp, col= 1:3, label = TRUE)

# same in ggplot2
NMDS = data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2])
ggplot(data = NMDS, aes(MDS1, MDS2)) + 
  geom_point(aes(data = meta_table, color = meta_table$Host_spp))

# First, make NMDS data frame with group column.

NMDS = data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2],group=meta_table$Host_spp)

# Next, save result of function ordiellipse() as some object.

ord<-ordiellipse(sol, meta_table$Host_spp, display = "sites", 
                 kind = "se", conf = 0.95, label = T)

#Data frame df_ell contains values to show ellipses. It is calculated again with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and now it uses arguments stored in ord object - cov, center and scale of each level.
#Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#Generate ellipse points

df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}
# Plotting is done the same way as in previous example. As for the calculating of coordinates for elipses object of ordiellipse() is used, this solution will work with different parameters you provide for this function.

ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = group)) +
  geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)

#UPDATE - solution that works in both cases
#First, make NMDS data frame with group column.
# $Host_spp

NMDS = data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2],group=meta_table$Host_spp)
#Next, save result of function ordiellipse() as some object.

ord<-ordiellipse(sol, meta_table$Host_spp, display = "sites", 
                 kind = "se", conf = 0.95, label = T)
ordispider(benthic.mds, meta_table$Host_spp, col= 1:3, label = TRUE)
#Data frame df_ell contains values to show ellipses. It is calculated again with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and now it uses arguments stored in ord object - cov, center and scale of each level.
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}
#Plotting is done the same way as in previous example. As for the calculating of coordinates for elipses object of ordiellipse() is used, this solution will work with different parameters you provide for this function.

ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = group)) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=1, linetype=2)

# $Host_size
NMDS = data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2],group=meta_table$Host_size)
#Next, save result of function ordiellipse() as some object.

ord<-ordiellipse(sol, meta_table$Host_size, display = "sites", 
                 kind = "se", conf = 0.95, label = T)
ordispider(sol, meta_table$Host_size, col= 1:3, label = TRUE)
#Data frame df_ell contains values to show ellipses. It is calculated again with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and now it uses arguments stored in ord object - cov, center and scale of each level.
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}
#Plotting is done the same way as in previous example. As for the calculating of coordinates for elipses object of ordiellipse() is used, this solution will work with different parameters you provide for this function.

ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = group)) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=1, linetype=2)
