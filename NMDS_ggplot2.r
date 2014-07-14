install.packages("vegan")
install.packages("colorspace")
install.packages("ggplot2")
install.packages("akima")
install.packages("Cairo")


library(vegan)
library(ggplot2)
library(ellipse)
library(grid)
library(akima)
library(reshape2)
library(Cairo)


# Load data file
x <-read.csv("C:/Users/Daniel/Dropbox/Projects/AZO/Analysis/AZO_otus_NMDS_worm_only.csv", header=TRUE)
x.env <-read.csv("C:/Users/Daniel/Dropbox/Projects/AZO/Analysis/AZO_map_worm_only.csv", header=TRUE)

#Transpose so its one rows:sample, column:OTU
x_matrix <- t(x)


#Calculate distance
x.dis <- vegdist(x_matrix)

#View stresses
x.mds0 <-monoMDS(x.dis)
stressplot(x.mds0, x.dis)


#Make MDS
x.mds <- metaMDS(x_matrix, trace = FALSE)

plot(x.mds, type = "p")		##t = labeled, p = points
fac1<-x.env$OriginDest			## choose factor for grouping
fac2<-x.env$LineA
ordispider(x.mds,group=fac1, label=FALSE, lwd=1)
ordiellipse(x.mds,group=fac1, label=TRUE, lwd=5)

ef <- envfit(x.mds ~ Al_soil + P_soil + S_soil + Ca_soil + Ti_soil + V_soil + Cr_soil + Mn_soil + Fe_soil + Co_soil + Ni_soil + Cu_soil + Zn_soil + As_soil + Se_soil + Mo_soil + Cd_soil + Sb_soil + Ba_soil + Hg_soil + Pb_soil + pH_soil + moi_soil + OC_soil, na.rm = TRUE, x.env)
plot(ef)
###############################################
#####             ggplot             ##########
###############################################

#Extract NMDS coordinates
sites <- scores(x.mds, display = "sites")

#Merge coordinates with env data
row.names(x.env) <-x.env$Sample
all_merge <- merge(sites, x.env, by=0, all=TRUE)
row.names(all_merge) <-all_merge$SampleID

all_merge

##interpolation##
d1 <- with(all_merge, interp(x = NMDS1, y = NMDS2, z = shannon, xo = seq(min(NMDS1), max(NMDS1), length = 360), yo = seq(min(NMDS2), max(NMDS2), length = 360)))
d2 <- melt(d1$z, na.rm = TRUE)
names(d2) <- c("x", "y", "shannon")
d2$NMDS1 <- d1$x[d2$x]
d2$NMDS2 <- d1$y[d2$y]

#### Environmental fits ####
ef_ori <- envfit(x.mds ~ ori_Altitude + ori_Temp , na.rm = TRUE, all_merge)
ef_dest <- envfit(x.mds ~ dest_Altitude + dest_Temp , na.rm = TRUE, all_merge)
ef_worm <- envfit(x.mds ~ worm_Mg + worm_Al, na.rm = TRUE, all_merge)

ef_all <- envfit(x.mds ~ ori_Altitude, na.rm = TRUE, all_merge)

ef_ori.df<-as.data.frame(ef_ori$vectors$arrows*sqrt(ef_ori$vectors$r))
ef_ori.df$species<-rownames(ef_ori.df)
ef_ori.df

ef_dest.df<-as.data.frame(ef_dest$vectors$arrows*sqrt(ef_dest$vectors$r))
ef_dest.df$species<-rownames(ef_dest.df)
ef_dest.df

ef_worm.df<-as.data.frame(ef_worm$vectors$arrows*sqrt(ef_worm$vectors$r))
ef_worm.df$species<-rownames(ef_worm.df)
ef_worm.df

ef_all.df<-as.data.frame(ef_all$vectors$arrows*sqrt(ef_all$vectors$r))
ef_all.df$species<-rownames(ef_all.df)
ef_all.df




##plots##
#### COMBINATION PLOT ####
 p <- ggplot() + theme_classic()  ## make blank plot as we'll be combining different data sets

pdf("/Users/Dan/Dropbox/Projects/DGC/figures/improved_NMDS.pdf", width = 28, height = 18)
CairoPNG(filename = "/Users/Dan/Dropbox/Projects/DGC/figures/improved_NMDS.png", width = 2048, height = 1024)
p + 
   stat_contour(data = d2, alpha=0.5, aes(x = NMDS1, y = NMDS2, fill = shannon, z = shannon)) +
   geom_raster(data = d2, alpha=0.5, aes(x = NMDS1, y = NMDS2, fill = exp(shannon), z = shannon)) + scale_fill_gradient(low="white", high="black") +
   geom_point(data=all_merge, size = 6, aes(NMDS1,NMDS2, shape = factor(Destination), colour = factor(Origin))) +

   theme(text = element_text(size = 10))

dev.off()
###################

##Overlay
   geom_segment(data=ef_all.df, aes(x=0,xend=NMDS1,y=0,yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")),colour="blue", size = 1, inherit_aes=FALSE) +
   geom_text(data=ef_all.df, aes(NMDS1,NMDS2,label=species),size=5, colour="black") + 

   geom_segment(data=ef_dest.df, aes(x=0,xend=NMDS1,y=0,yend=NMDS2), arrow = arrow(length = unit(0.5, "cm")),colour="red", size = 1, inherit_aes=FALSE) +
   geom_segment(data=ef_worm.df, aes(x=0,xend=NMDS1,y=0,yend=NMDS2), arrow = arrow(length = unit(0.5, "cm")),colour="green", size = 1, inherit_aes=FALSE) +
   geom_segment(data=ef_ori.df, aes(x=0,xend=NMDS1,y=0,yend=NMDS2), arrow = arrow(length = unit(0.5, "cm")),colour="blue", size = 1, inherit_aes=FALSE) +

   geom_text(data=ef_dest.df, aes(NMDS1,NMDS2,label=species),size=3, colour="black") + 
   geom_text(data=ef_worm.df, aes(NMDS1,NMDS2,label=species),size=3, colour="black") + 
   geom_text(data=ef_ori.df, aes(NMDS1,NMDS2,label=species),size=3, colour="black") + 


###basic plots
 p + geom_point(size = 4, alpha=.8, aes(colour = factor(Site)))
 p + geom_point(size = 4, alpha=.8, aes(colour = factor(Site))) + geom_path(data=df_ell, aes(x=x, y=y, colour=group), size=1, linetype=1)
 p + geom_point(size = 4, alpha=.8, aes(colour = factor(Condition), shape = factor(Condition))) + geom_path(data=df_ell, aes(x=x, y=y, colour=group), size=1, linetype=1)



p + geom_point(size = 4, alpha=.8, aes(colour = factor(Site), shape = factor(Astype))) + 
   geom_segment(data=ef.df,aes(x=0,xend=NMDS1,y=0,yend=NMDS2), arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE) +
   geom_text(data=ef.df,aes(x=NMDS1,y=NMDS2,label=species),size=5)
      

## With ellipse##

##Ellipses##
df_ell <- data.frame()
df_ell

for(g in levels(all_merge$Site)){
df_ell <- rbind(df_ell, cbind(as.data.frame(with(all_merge[all_merge$Site==g,], ellipse(cor(NMDS1,NMDS2), 
        scale=c(sd(NMDS1),sd(NMDS2)), 
        centre=c(mean(NMDS1),mean(NMDS2))))),group=g))}


p + geom_point(size = 4, alpha=.8, aes(colour = factor(Site), shape = factor(Site))) + 
   geom_path(data=df_ell, aes(x=x, y=y, colour=group), size=1, linetype=1) + 
   geom_segment(data=ef.df,aes(x=0,xend=NMDS1,y=0,yend=NMDS2), arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE) + 
   geom_text(data=ef.df,aes(x=NMDS1,y=NMDS2,label=Species),size=5)

