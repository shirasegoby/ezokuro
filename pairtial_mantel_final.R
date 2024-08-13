getwd()
setwd("/Users/shotarohirase/Desktop/エゾクロ解析re/IBDIBE_plot/")

#par(mfrow=c(1,1)) 
# Load package
#install.packages("marmap")
#install.packages("tcR")
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
#install.packages("tcr")
#install.packages("immunarch") 
#library(immunarch)

Lon.range = c(128, 147)
Lat.range = c(30, 46)
bat <- getNOAA.bathy(Lon.range[1],Lon.range[2],Lat.range[1],Lat.range[2],res=5)

#Load land and bathymetry data
trans1 <- trans.mat(bat)

#******calculate geographic distance**********
# Create vectors of latitude and longitude
lon <- c(140.710137,140.995927,141.578335,142.071,142.002381,141.831514,141.545721,141.0336,140.631186,140.0368365,139.615667,138.94444,136.926077,136.744246,134.45030536,131.132496,141.140754,141.50257379990148,140.101929,140.238494,139.6394821,139.829594,139.549427,139.58706,139.258123,139.425016,138.322569,136.107972,133.34107,129.95619,128.865569)
lat <- c(42.520818,41.52102,40.544865,39.731328,39.266691,38.998739,38.3525,36.968269,36.377239,34.974629,35.128866,34.647859,34.442752,34.169183,33.628055,31.438253,45.249765,44.00240811743948,41.420012,41.132959,39.9795329,39.2325526,39.194308,38.63859,38.46598,38.188655,37.964991,36.231662,36.156322,33.548338,32.711916)
#points(lon, lat, pch = 21, bg = "red2", cex = 0.8)

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

# PO vs SJ
pdf("IBD_FST_GRAS.pdf", width = 10, height = 10)
par(mar = c(10, 5, 10, 5))  # マージンを調整
#地理的距離とFSTのプロット
plot(dist1,fst.dist,ylim=c(0.00,0.30),xlim=c(0,2000))
#abline(lm(dist1~fst.dist),col="red",lty=2)
#PO vs SJ
par(new = TRUE)
plot(as.matrix(dist1)[17:31,1:16],as.matrix(fst.dist)[17:31,1:16],col = "blue",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))
# PO vs PO
par(new = TRUE)
plot(as.matrix(dist1)[1:16,1:16],as.matrix(fst.dist)[1:16,1:16],col = "orange",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))
# SJ vs SJ
par(new = TRUE)
plot(as.matrix(dist1)[17:31,17:31],as.matrix(fst.dist)[17:31,17:31],col = "green",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))
dev.off()

# ezo v kuro
pdf("IBD_FST_GRAS_ezokuro.pdf", width = 10, height = 10)
par(mar = c(10, 5, 10, 5))  
#geo distance vs FST
plot(dist1,fst.dist,ylim=c(0.00,0.30),xlim=c(0,2000))
#abline(lm(dist1~fst.dist),col="red",lty=2)
#ezo vs kuro
par(new = TRUE)
plot(as.matrix(dist1)[1:9,10:16],as.matrix(fst.dist)[1:9,10:16],col = "blue",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))
par(new = TRUE)
plot(as.matrix(dist1)[1:9,20:31],as.matrix(fst.dist)[1:9,20:31],col = "blue",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))
par(new = TRUE)
plot(as.matrix(dist1)[10:16,17:19],as.matrix(fst.dist)[10:16,17:19],col = "blue",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))
par(new = TRUE)
plot(as.matrix(dist1)[17:19,20:31],as.matrix(fst.dist)[17:19,20:31],col = "blue",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))

#ezo
par(new = TRUE)
plot(as.matrix(dist1)[1:9,1:9],as.matrix(fst.dist)[1:9,1:9],col = "green",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))
par(new = TRUE)
plot(as.matrix(dist1)[1:9,17:19],as.matrix(fst.dist)[1:9,17:19],col = "green",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))
par(new = TRUE)
plot(as.matrix(dist1)[17:19,17:19],as.matrix(fst.dist)[17:19,17:19],col = "green",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))

#kuro
par(new = TRUE)
plot(as.matrix(dist1)[10:16,10:16],as.matrix(fst.dist)[10:16,10:16],col = "orange",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))
par(new = TRUE)
plot(as.matrix(dist1)[10:16,20:31],as.matrix(fst.dist)[10:16,20:31],col = "orange",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))
par(new = TRUE)
plot(as.matrix(dist1)[20:31,20:31],as.matrix(fst.dist)[20:31,20:31],col = "orange",pch=16,ylim=c(0.00,0.30),xlim=c(0,2000))

dev.off()


# data convert
#ezokuro
data1 <- as.matrix(dist1)[1:7,8:16]
data2 <- as.matrix(fst.dist)[1:7,8:16]
ezokuro1 <- data.frame(geo = as.vector(data1), fst = as.vector(data2))
data1 <- as.matrix(dist1)[1:7,20:31]
data2 <- as.matrix(fst.dist)[1:7,20:31]
ezokuro2 <- data.frame(geo = as.vector(data1), fst = as.vector(data2))
data1 <- as.matrix(dist1)[8:16,17:19]
data2 <- as.matrix(fst.dist)[8:16,17:19]
ezokuro3 <- data.frame(geo = as.vector(data1), fst = as.vector(data2))
data1 <- as.matrix(dist1)[17:19,20:31]
data2 <- as.matrix(fst.dist)[17:19,20:31]
ezokuro4 <- data.frame(geo = as.vector(data1), fst = as.vector(data2))

ezokuro <- rbind(ezokuro1,ezokuro2,ezokuro3,ezokuro4)
ezokuro_model <- lmodel2(formula =ezokuro$fst~ezokuro$geo)

par(new = TRUE)
plot(ezokuro_model,"SMA",xlim=c(0,2000),ylim=c(0,0.30), xlab="", ylab="",pch = 16,col=c("red","orange","orange"))
line.cis(ezokuro$fst,ezokuro$geo)

#ezo
data1 <- as.matrix(dist1)[1:7,1:7]
data2 <- as.matrix(fst.dist)[1:7,1:7]
ezo1 <- data.frame(geo = as.vector(data1), fst = as.vector(data2))
data1 <- as.matrix(dist1)[1:7,17:19]
data2 <- as.matrix(fst.dist)[1:7,17:19]
ezo2 <- data.frame(geo = as.vector(data1), fst = as.vector(data2))
data1 <- as.matrix(dist1)[17:19,17:19]
data2 <- as.matrix(fst.dist)[17:19,17:19]
ezo3 <- data.frame(geo = as.vector(data1), fst = as.vector(data2))

ezo <- rbind(ezo1,ezo2,ezo3)
ezo_model <- lmodel2(formula =ezo$fst~ezo$geo)

par(new = TRUE)
plot(ezo_model,"SMA",xlim=c(0,2000),ylim=c(0,0.30), xlab="", ylab="",pch = 16,col=c("yellow","yellow","yellow"))
line.cis(ezo$fst,ezo$geo)

#kuro
data1 <- as.matrix(dist1)[8:16,8:16]
data2 <- as.matrix(fst.dist)[8:16,8:16]
kuro1 <- data.frame(geo = as.vector(data1), fst = as.vector(data2))

data1 <- as.matrix(dist1)[8:16,20:31]
data2 <- as.matrix(fst.dist)[8:16,20:31]
kuro2 <- data.frame(geo = as.vector(data1), fst = as.vector(data2))

data1 <- as.matrix(dist1)[20:31,20:31]
data2 <- as.matrix(fst.dist)[20:31,20:31]
kuro3 <- data.frame(geo = as.vector(data1), fst = as.vector(data2))

kuro <- rbind(kuro1,kuro2,kuro3)
kuro_model <- lmodel2(formula =kuro$fst~kuro$geo)

par(new = TRUE)
plot(kuro_model,"SMA",xlim=c(0,2000),ylim=c(0,0.30), xlab="", ylab="",pch = 16,col=c("red","orange","orange"))
line.cis(kuro$fst,kuro$geo)

#slope
ezokuro_model$regression.results[3, ] %>% kable()
kuro_model$regression.results[3, ] %>% kable()
ezo_model$regression.results[3, ] %>% kable()

#group列の追加
ezokuro$group <- rep("ezokuro", nrow(ezokuro))
kuro$group <- rep("kuro", nrow(kuro))
ezo$group <- rep("ezo", nrow(ezo))

#全部のデータを合わせる
all_data <- rbind(ezokuro,kuro,ezo)


mod3 <- aov (fst ~ geo*group, data = all_data)
summary(mod3)


ezo_model <- lm(fst ~ geo, data = ezo)
kuro_model <- lm(fst ~ geo, data = kuro)
ezokuro_model <- lm(fst ~ geo, data = ezokuro)

model_comparison <- anova(kuro_model,ezo_model,ezokuro_model)

#*********calculate env data****************
# env data
all <- read.table("env.txt",header=T)
rownames(all)
#PCA解析
all_pca <- prcomp(all, scale = TRUE)
#biplot(x=all_pca)
#PC1の貢献度を%を確認する
summary(all_pca)
#PC1    PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10   PC11    PC12
#Standard deviation     4.4290 2.0784 1.53880 1.31258 1.07263 0.93075 0.91138 0.79109 0.66272 0.54404 0.4259 0.41401
#Proportion of Variance 0.5944 0.1309 0.07175 0.05221 0.03486 0.02625 0.02517 0.01896 0.01331 0.00897 0.0055 0.00519
#Cumulative Proportion  0.5944 0.7253 0.79708 0.84928 0.88415 0.91040 0.93557 0.95454 0.96784 0.97681 0.9823 0.98750
#PC13    PC14    PC15    PC16    PC17    PC18    PC19    PC20    PC21    PC22    PC23    PC24
#Standard deviation     0.33923 0.29564 0.26458 0.19971 0.18733 0.13811 0.11712 0.10111 0.07481 0.06530 0.06466 0.05924
#Proportion of Variance 0.00349 0.00265 0.00212 0.00121 0.00106 0.00058 0.00042 0.00031 0.00017 0.00013 0.00013 0.00011
#Cumulative Proportion  0.99099 0.99364 0.99576 0.99697 0.99803 0.99861 0.99903 0.99934 0.99951 0.99964 0.99976 0.99987
#PC25    PC26    PC27    PC28    PC29    PC30      PC31
#Standard deviation     0.04553 0.02698 0.02592 0.02165 0.01526 0.01231 1.734e-15
#Proportion of Variance 0.00006 0.00002 0.00002 0.00001 0.00001 0.00000 0.000e+00
#Cumulative Proportion  0.99993 0.99995 0.99997 0.99999 1.00000 1.00000 1.000e+00


##環境のPCAプロットを作成する
##PC1が環境の南北差を反映しているか確認する
p1.pc <- as.data.frame(all_pca$x)
p1.pc$pop <- c("PO_1","PO_2","PO_3","PO_4","PO_5","PO_6","PO_7","PO_8","PO_9","PO_10","PO_11","PO_12","PO_13","PO_14","PO_15","PO_16","SJ_1","SJ_2","SJ_3","SJ_4","SJ_5","SJ_6","SJ_7","SJ_8","SJ_9","SJ_10","SJ_11","SJ_12","SJ_13","SJ_14","SJ_15")
#indexを固定
p1.pc$pop <- factor(p1.pc$pop, levels = c("PO_1","PO_2","PO_3","PO_4","PO_5","PO_6","PO_7","PO_8","PO_9","PO_10","PO_11","PO_12","PO_13","PO_14","PO_15","PO_16","SJ_1","SJ_2","SJ_3","SJ_4","SJ_5","SJ_6","SJ_7","SJ_8","SJ_9","SJ_10","SJ_11","SJ_12","SJ_13","SJ_14","SJ_15")) #これ！！！！！！！！！！！！！！！！！！！
rownames(p1.pc) <- c("PO_1","PO_2","PO_3","PO_4","PO_5","PO_6","PO_7","PO_8","PO_9","PO_10","PO_11","PO_12","PO_13","PO_14","PO_15","PO_16","SJ_1","SJ_2","SJ_3","SJ_4","SJ_5","SJ_6","SJ_7","SJ_8","SJ_9","SJ_10","SJ_11","SJ_12","SJ_13","SJ_14","SJ_15")

#color setting
colors <- colorRampPalette(c("light blue","blue", "purple","green", "orange","red"))

#PA vs SJ
col_list <- c(colors(16),colors(15))

#PC1 and PC2 plot
#個体ラベルあり
g <- ggplot() +
  #geom_point(data=p1.pc, aes(x=PC1, y=PC2, shape=Species, colour=Species), size=4)+
  geom_point(data=p1.pc, aes(x=PC1, y=PC2,color=pop), size=4)+
  scale_color_manual(values = col_list)+
  ggrepel::geom_label_repel(data =p1.pc, 
                            mapping = aes(x=PC1, y=PC2, label = rownames(p1.pc)))
g

#PC1 score
PC1score <- all_pca$x[,1] 
#PC2 score
PC2score <- all_pca$x[,2] 

data <- data.frame(p1.pc$pop,PC1score)
data2 <- data.frame(p1.pc$pop,PC2score)
data

plot.new()
par(mar = c(10, 5, 5, 5))

#背景を透明に
par(bg = "transparent")

x <- c(1:31)

pdf("PC1_env.pdf", width = 10, height = 8)
beeswarm(data = data, PC1score ~ p1.pc$pop, las = 2, col = col_list, pch = 16,cex=2.5　, xlab = "", ylab="PC1 score")  # プロットを追加
abline(v = x, col = "grey", lty = 3)  # ablineをプロットの下に配置
par(bg = "transparent")
par(new=T)
beeswarm(data = data, PC1score ~ p1.pc$pop, las = 2, col = col_list, pch = 16,cex=2.5, xlab = "")
dev.off()


pdf("PC2_env.pdf", width = 13, height = 5)
beeswarm(data = data, PC2score ~ p1.pc$pop, las = 2, col = col_list, pch = 16, xlab = "", ylab="PC1 score")  # プロットを追加
abline(v = x, col = "grey", lty = 3)  # ablineをプロットの下に配置
par(bg = "transparent")
par(new=T)
beeswarm(data = data, PC2score ~ p1.pc$pop, las = 2, col = col_list, pch = 16, xlab = "")
dev.off()

##IBE plot
#make distance matrix
all_env_dist <- dist(PC1score, method = "euclidean")
#env distance
ibe <- vegan::mantel(fst.dist,all_env_dist) #veganのMatel
sink ("env_Mantel_result.txt")
ibe
sink ()

#env distance vs FST plot
pdf("IBE_FST_GRAS.pdf", width = 10, height = 10)
par(mar = c(10, 5, 10, 5)) 

plot(all_env_dist,fst.dist,ylim=c(0.00,0.30),xlim=c(0,16))
#abline(lm(PO_env_dist~PO_fst.dist),col="red",lty=2)
par(new = TRUE)
plot(as.matrix(all_env_dist)[17:31,1:16],as.matrix(fst.dist)[17:31,1:16],col = "blue",pch=16,ylim=c(0.00,0.30),xlim=c(0,16))

# PO vs PO
par(new = TRUE)
plot(as.matrix(all_env_dist)[1:16,1:16],as.matrix(fst.dist)[1:16,1:16],col = "orange",pch=16,ylim=c(0.00,0.30),xlim=c(0,16))

# SJ vs SJ
par(new = TRUE)
plot(as.matrix(all_env_dist)[17:31,17:31],as.matrix(fst.dist)[17:31,17:31],col = "green",pch=16,ylim=c(0.00,0.30),xlim=c(0,16))
dev.off()

#***************partial Mantel test********************
partial <- vegan::mantel.partial(fst.dist,all_env_dist,dist1) #veganのMatel
sink ("partial_Mantel_result.txt")
partial
sink ()






