

library(pcadapt)
library("poppr")
library("PopGenReport")
library("adegenet")
library("genepopedit")
library(vioplot)
library(beeswarm)

#＃https://cran.r-project.org/web/packages/pcadapt/vignettes/pcadapt.html

getwd()
setwd("./pcadapt")

data <- read.pcadapt("ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped",type="ped")

ind <- read.table("ezokuro_filtered.recode.dDocent.recode.rename2.sorted.ped",header=F)

poplist.names <- ind$V1
poplist.names

#PC number
K <- 10
# run
x <- pcadapt(data,K=K)
summary(x)

#contribution
(x$singular.values)^2

#popの名前を取得
ind2 <- read.table("ind_category",header=F)

#V1が、popname
names<- ind2$V1

PC1score <- x$scores[,1]
PC2score <- x$scores[,2]

#カラーを指定
colors <- colorRampPalette(c("light blue","blue", "purple","green", "orange","red"))

#太平洋と日本海それぞれ
col_list <- c(colors(16),colors(15))

plot(x,option="scores",pop=ind2$V2,col=col_list)
plot(x,option="scores",i=3,j=4,pop=ind2$V2,col=col_list)

# data modify
data <- data.frame(ind2$V1,PC1score)
data2 <- data.frame(ind2$V1,PC2score)

plot.new()
par(mar = c(10, 5, 5, 5))

par(bg = "transparent")

x <- c(1:31)

pdf("PC1.pdf", width = 13, height = 5)
beeswarm(data = data, PC1score ~ ind2$V1, las = 2, col = col_list, pch = 16, xlab = "", ylab="PC1 score") 
abline(v = x, col = "grey", lty = 3) 
par(bg = "transparent")
par(new=T)
beeswarm(data = data, PC1score ~ ind2$V1, las = 2, col = col_list, pch = 16, xlab = "")
dev.off()


plot.new()
par(mar = c(10, 5, 5, 5))

par(bg = "transparent")

x <- c(1:31)

pdf("PC2.pdf", width = 13, height = 5)
beeswarm(data = data2, PC2score ~ ind2$V1, las = 2, col = col_list, pch = 16, xlab = "", ylab="PC2 score")  # プロットを追加
abline(v = x, col = "grey", lty = 3)  # ablineをプロットの下に配置
par(bg = "transparent")
par(new=T)
beeswarm(data = data2, PC2score ~ ind2$V1, las = 2, col = col_list, pch = 16, xlab = "")
dev.off()


#########outlier SNP detection######
locus <- read.table("ezokuro_filtered.recode.dDocent.recode.rename2.sorted.map",header=F)
#locus$V1
#locus$V2
#locus$V4

x <- pcadapt(data,K=1)
summary(x)
plot(x,option="manhattan")
#x$pvalues

#FDR cutoff
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
p_q_value <-cbind(as.character(locus$V1),as.character(locus$V2),as.character(locus$V4),x$pvalues,qval)

#threthold
alpha <- 0.05

all <- which(qval<1)

snp_pc <- get.pc(x,all)
out <- cbind(na.omit(p_q_value),snp_pc$PC)

out <-data.frame(out)

out2 <- subset(out,as.numeric(as.character(out$qval))<0.05 & out$V6==1)
write.table(out2,"pcadapt_PC1_outlier.txt",quote=F,col.names=F,row.names=F)


