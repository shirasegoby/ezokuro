getwd()
setwd("./IBDIBE_plot/")

install.packages("smatr")
library(lmodel2)
library(gdistance)
library(adegenet)
library(PopGenReport)
library(hierfstat)
library(marmap)
library(vegan)
library(ggplot2)
library(beeswarm)
library(smatr)

Lon.range = c(128, 147)
Lat.range = c(30, 46)
bat <- getNOAA.bathy(Lon.range[1],Lon.range[2],Lat.range[1],Lat.range[2],res=5)

#Load land and bathymetry data
trans1 <- trans.mat(bat)

#******calculate geographic distance**********
# Create vectors of latitude and longitude
lon <- c(140.710137,140.995927,141.578335,142.071,142.002381,141.831514,141.545721,141.0336,140.631186,140.0368365,139.615667,138.94444,136.926077,136.744246,134.45030536,131.132496,141.140754,141.50257379990148,140.101929,140.238494,139.6394821,139.829594,139.549427,139.58706,139.258123,139.425016,138.322569,136.107972,133.34107,129.95619,128.865569)
lat <- c(42.520818,41.52102,40.544865,39.731328,39.266691,38.998739,38.3525,36.968269,36.377239,34.974629,35.128866,34.647859,34.442752,34.169183,33.628055,31.438253,45.249765,44.00240811743948,41.420012,41.132959,39.9795329,39.2325526,39.194308,38.63859,38.46598,38.188655,37.964991,36.231662,36.156322,33.548338,32.711916)

sites <- cbind(lon,lat)
rownames(sites) <- 1:31
#sites
#lc.dist: Computes least cost distances between two or more locations
dist1 <- marmap::lc.dist(trans1, sites, res = "dist")

# Create nice looking color palettes
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

# map plot
plot(bat, image = TRUE, land = TRUE, lwd = 0.1, bpal = list(c(0, max(bat), greys), c(min(bat), 0, blues)))
plot(bat, lwd = 0.8, deep = 0, shallow = 0, step = 0, add = TRUE) # highlight coastline

points(lon, lat, pch = 21, bg = "red2", cex = 0.8)

#*****FST value****************************
genedis <- read.csv("pairFST_genodive.csv",header=T)
rownames(genedis) <- genedis$Statistic
genedis <- genedis[,-1]

#dist object
fst.dist <- as.dist(genedis)

#FST/1-FST
fst.dist2 <- (fst.dist/(1-fst.dist))
#write.table(genedis2,"SJ_WC_FST.txt")

#************Mantel test*****************
#geographic distance
ibd <-  vegan::mantel(fst.dist,dist1) #veganのMatel
sink ("dis_Mantel_result.txt")
ibd
sink ()

#*********calculate env data****************
# env data
all <- read.table("env.txt",header=T)
rownames(all)
#PCA
all_pca <- prcomp(all, scale = TRUE)
summary(all_pca)

#PC1 score
PC1score <- all_pca$x[,1] 
#PC2 score
PC2score <- all_pca$x[,2] 

data <- data.frame(p1.pc$pop,PC1score)
data2 <- data.frame(p1.pc$pop,PC2score)

#make distance matrix
all_env_dist <- dist(PC1score, method = "euclidean")
#env distance
ibe <- vegan::mantel(fst.dist,all_env_dist) 
sink ("env_Mantel_result.txt")
ibe
sink ()

#***************partial Mantel test********************
partial <- vegan::mantel.partial(fst.dist,all_env_dist,dist1) #veganのMatel
sink ("partial_Mantel_result.txt")
partial
sink ()






