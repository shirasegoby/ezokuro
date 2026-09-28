setwd("./outflank/")

library(devtools)
library(LEA)
library("OutFLANK")
library("poppr")
library("PopGenReport")
library("adegenet")
library("genepopedit")

if (!("devtools" %in% installed.packages())){install.packages(devtools)}
library(devtools)

if (!("qvalue" %in% installed.packages())){TODO}
if (!("vcfR" %in% installed.packages())){install.packages("vcfR")} 

data <- read.table("ezokuro_filtered.recode.dDocent.recode.rename2.sorted.OUTflank2.raw",header=F,skip=1)

data2 <- data[, -which (colnames(data) %in% c("V1", "V2", "V3","V4","V5","V6"))]

write.table(data2,"mod.ped",quote=F,col.names=F,row.names=F)

system("sed -e s/NA/9/g mod.ped > mod2.ped")

SNPmat <- read.table("mod2.ped",header=F)

locusNames <- read.table("ezokuro_filtered.recode.dDocent.recode.rename2.sorted.map",header=F)
length(locusNames$V1)
popNames <- c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5),rep(7,5),rep(8,5),rep(9,5),rep(10,5),rep(11,5),
              rep(12,5),rep(13,5),rep(14,5),rep(15,5),rep(16,5),rep(17,5),rep(18,5),rep(19,5),rep(20,5),rep(21,5),rep(22,5),
              rep(23,5),rep(24,5),rep(25,5),rep(26,4),rep(27,5),rep(28,5),rep(29,5),rep(30,5),rep(31,5))
length(popNames)

abalone <- MakeDiploidFSTMat(SNPmat,locusNames$V2,popNames)

#relationship between He and FST
plot(abalone$He, abalone$FSTNoCorr, pch=10, col="grey")

hist(abalone$FSTNoCorr, breaks=seq(0,1, by=0.01))

# Removing low Heterozygosity variants results in a more chi-square looking FST distribution
hist(abalone$FSTNoCorr[abalone$He>0.1], breaks=seq(0,1, by=0.01))
#He > 0.1 is good

###running outflank
# increasing the Right Trim Fraction doesn't help, but increasing the left trim fraction enables a better fit
#following http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html
outflank_result <- OutFLANK(FstDataFrame=abalone, NumberOfSamples=31,Hmin=0.1,
                            qthreshold=0.05,RightTrimFraction = 0.06, LeftTrimFraction = 0.35)

hist(outflank_result$results$pvaluesRightTail)

OutFLANKResultsPlotter(outflank_result, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

OutFLANKResultsPlotter(outflank_result, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

OutFLANKResultsPlotter(outflank_result, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.4, titletext = NULL)

#write result once
write.csv(outflank_result,"outflank_result.csv")

true_site <- subset(outflank_result$results,outflank_result$results$OutlierFlag=="TRUE")

write.csv(true_site,"outflank_result_outlier.csv")


#outlier number
dim(true_site)[1]

#all locus from map file
#subset outlier
#subset(locusNames, V2 %in% true_site$LocusName)
write.table(subset(locusNames, V2 %in% true_site$LocusName),"outflank_outlier.map",quote=F,col.names=F,row.names=F)

