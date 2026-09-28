# Isolation-by-distance (IBD) and isolation-by-environment (IBE) analyses
#
# Required R packages:
#   lmodel2, gdistance, adegenet, PopGenReport, hierfstat,
#   marmap, vegan, ggplot2, ggrepel, beeswarm, smatr
#
# Required input files:
#   pairFST_genodive.csv
#   env_md.txt
#
# The script assumes that input files are located in the working directory.


########################################
# 1. Load packages
########################################

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


########################################
# 2. Bathymetric data
########################################

lon_range <- c(128, 147)
lat_range <- c(30, 46)

bat <- getNOAA.bathy(
    lon_range[1],
    lon_range[2],
    lat_range[1],
    lat_range[2],
    res = 5
)

# Convert bathymetric data to a transition matrix.
trans1 <- trans.mat(bat)


########################################
# 3. Geographic distances among sites
########################################

# Sampling-site coordinates.
lon <- c(
    140.71014, 140.99593, 141.57834, 142.07100, 142.00238,
    141.83151, 141.54572, 141.03360, 140.63119, 140.03684,
    139.61567, 138.94444, 136.92608, 136.74425, 134.45031,
    131.13250, 141.14075, 141.50257, 140.10193, 140.23849,
    139.63948, 139.82959, 139.55491, 139.58706, 139.25812,
    138.32257, 139.42502, 136.10797, 133.34107, 129.95619,
    128.86557
)

lat <- c(
    42.52082, 41.52102, 40.54487, 39.73133, 39.26669,
    38.99874, 38.35250, 36.96827, 36.37724, 34.97463,
    35.12887, 34.64786, 34.44275, 34.16918, 33.62806,
    31.43825, 45.24977, 44.00241, 41.42001, 41.13296,
    39.97953, 39.23255, 39.18958, 38.63859, 38.46598,
    37.96499, 38.18866, 36.23166, 36.15632, 33.54834,
    32.71192
)

sites <- cbind(lon, lat)
rownames(sites) <- 1:31

# Least-cost geographic distances following marine pathways.
dist1 <- marmap::lc.dist(trans1, sites, res = "dist")

# Optional:
# write.csv(as.matrix(dist1), "geographic_distance.csv")


########################################
# 4. Plot sampling locations
########################################

blues <- c(
    "lightsteelblue4",
    "lightsteelblue3",
    "lightsteelblue2",
    "lightsteelblue1"
)

greys <- c(grey(0.6), grey(0.93), grey(0.99))

plot(
    bat,
    image = TRUE,
    land = TRUE,
    lwd = 0.1,
    bpal = list(
        c(0, max(bat), greys),
        c(min(bat), 0, blues)
    )
)

plot(
    bat,
    lwd = 0.8,
    deep = 0,
    shallow = 0,
    step = 0,
    add = TRUE
)

points(lon, lat, pch = 21, bg = "red2", cex = 0.8)


########################################
# 5. Pairwise FST
########################################

# Pairwise FST values were calculated using GenoDive.
genedis <- read.csv("pairFST_genodive.csv", header = TRUE)

rownames(genedis) <- genedis$Statistic
genedis <- genedis[, -1]

fst.dist <- as.dist(genedis)

# FST/(1-FST), retained from the original workflow.
fst.dist2 <- fst.dist / (1 - fst.dist)


########################################
# 6. Mantel test: isolation by distance
########################################

ibd <- vegan::mantel(fst.dist, dist1)

sink("dis_Mantel_result.txt")
print(ibd)
sink()


########################################
# 7. IBD plots
########################################

# PO vs. SJ comparison.
pdf("IBD_FST_GRAS.pdf", width = 10, height = 10)

par(mar = c(10, 5, 10, 5))

plot(
    dist1,
    fst.dist,
    ylim = c(0.00, 0.30),
    xlim = c(0, 2000)
)

# PO vs. SJ
par(new = TRUE)
plot(
    as.matrix(dist1)[17:31, 1:16],
    as.matrix(fst.dist)[17:31, 1:16],
    col = "blue",
    pch = 16,
    ylim = c(0.00, 0.30),
    xlim = c(0, 2000)
)

# PO vs. PO
par(new = TRUE)
plot(
    as.matrix(dist1)[1:16, 1:16],
    as.matrix(fst.dist)[1:16, 1:16],
    col = "orange",
    pch = 16,
    ylim = c(0.00, 0.30),
    xlim = c(0, 2000)
)

# SJ vs. SJ
par(new = TRUE)
plot(
    as.matrix(dist1)[17:31, 17:31],
    as.matrix(fst.dist)[17:31, 17:31],
    col = "green",
    pch = 16,
    ylim = c(0.00, 0.30),
    xlim = c(0, 2000)
)

dev.off()


# Ezo vs. Kuro comparison.
pdf("IBD_FST_GRAS_ezokuro.pdf", width = 10, height = 10)

par(mar = c(10, 5, 10, 5))

plot(
    dist1,
    fst.dist,
    ylim = c(0.00, 0.30),
    xlim = c(0, 2000)
)

# Ezo vs. Kuro
par(new = TRUE)
plot(as.matrix(dist1)[1:9, 10:16],
     as.matrix(fst.dist)[1:9, 10:16],
     col = "blue", pch = 16,
     ylim = c(0.00, 0.30), xlim = c(0, 2000))

par(new = TRUE)
plot(as.matrix(dist1)[1:9, 20:31],
     as.matrix(fst.dist)[1:9, 20:31],
     col = "blue", pch = 16,
     ylim = c(0.00, 0.30), xlim = c(0, 2000))

par(new = TRUE)
plot(as.matrix(dist1)[10:16, 17:19],
     as.matrix(fst.dist)[10:16, 17:19],
     col = "blue", pch = 16,
     ylim = c(0.00, 0.30), xlim = c(0, 2000))

par(new = TRUE)
plot(as.matrix(dist1)[17:19, 20:31],
     as.matrix(fst.dist)[17:19, 20:31],
     col = "blue", pch = 16,
     ylim = c(0.00, 0.30), xlim = c(0, 2000))

# Ezo
par(new = TRUE)
plot(as.matrix(dist1)[1:9, 1:9],
     as.matrix(fst.dist)[1:9, 1:9],
     col = "green", pch = 16,
     ylim = c(0.00, 0.30), xlim = c(0, 2000))

par(new = TRUE)
plot(as.matrix(dist1)[1:9, 17:19],
     as.matrix(fst.dist)[1:9, 17:19],
     col = "green", pch = 16,
     ylim = c(0.00, 0.30), xlim = c(0, 2000))

par(new = TRUE)
plot(as.matrix(dist1)[17:19, 17:19],
     as.matrix(fst.dist)[17:19, 17:19],
     col = "green", pch = 16,
     ylim = c(0.00, 0.30), xlim = c(0, 2000))

# Kuro
par(new = TRUE)
plot(as.matrix(dist1)[10:16, 10:16],
     as.matrix(fst.dist)[10:16, 10:16],
     col = "orange", pch = 16,
     ylim = c(0.00, 0.30), xlim = c(0, 2000))

par(new = TRUE)
plot(as.matrix(dist1)[10:16, 20:31],
     as.matrix(fst.dist)[10:16, 20:31],
     col = "orange", pch = 16,
     ylim = c(0.00, 0.30), xlim = c(0, 2000))

par(new = TRUE)
plot(as.matrix(dist1)[20:31, 20:31],
     as.matrix(fst.dist)[20:31, 20:31],
     col = "orange", pch = 16,
     ylim = c(0.00, 0.30), xlim = c(0, 2000))

dev.off()


# Ezo/Kuro and PO/SJ comparison.
pdf("IBD_FST_GRAS_ezokuroPOSJ.pdf", width = 10, height = 10)

par(mar = c(10, 5, 10, 5))

# pch = 24: PO vs. PO
# pch = 21: PO vs. SJ
# pch = 22: SJ vs. SJ

# Ezo PO vs. Kuro PO
plot(as.matrix(dist1)[1:9, 10:16],
     as.matrix(fst.dist)[1:9, 10:16],
     col = "grey", pch = 24,
     ylim = c(0.00, 0.30), xlim = c(0, 2000), cex = 1.5)

# Ezo PO vs. Kuro SJ
par(new = TRUE)
plot(as.matrix(dist1)[1:9, 20:31],
     as.matrix(fst.dist)[1:9, 20:31],
     col = "grey", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 2000), cex = 1.5)

# Kuro PO vs. Ezo SJ
par(new = TRUE)
plot(as.matrix(dist1)[10:16, 17:19],
     as.matrix(fst.dist)[10:16, 17:19],
     col = "grey", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 2000), cex = 1.5)

# Ezo SJ vs. Kuro SJ
par(new = TRUE)
plot(as.matrix(dist1)[17:19, 20:31],
     as.matrix(fst.dist)[17:19, 20:31],
     col = "grey", pch = 22,
     ylim = c(0.00, 0.30), xlim = c(0, 2000), cex = 1.5)

# Ezo PO vs. Ezo PO
par(new = TRUE)
plot(as.matrix(dist1)[1:9, 1:9],
     as.matrix(fst.dist)[1:9, 1:9],
     col = "blue", pch = 24,
     ylim = c(0.00, 0.30), xlim = c(0, 2000), cex = 1.5)

# Ezo PO vs. Ezo SJ
par(new = TRUE)
plot(as.matrix(dist1)[1:9, 17:19],
     as.matrix(fst.dist)[1:9, 17:19],
     col = "blue", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 2000), cex = 1.5)

# Ezo SJ vs. Ezo SJ
par(new = TRUE)
plot(as.matrix(dist1)[17:19, 17:19],
     as.matrix(fst.dist)[17:19, 17:19],
     col = "blue", pch = 22,
     ylim = c(0.00, 0.30), xlim = c(0, 2000), cex = 1.5)

# Kuro PO vs. Kuro PO
par(new = TRUE)
plot(as.matrix(dist1)[10:16, 10:16],
     as.matrix(fst.dist)[10:16, 10:16],
     col = "red", pch = 24,
     ylim = c(0.00, 0.30), xlim = c(0, 2000), cex = 1.5)

# Kuro PO vs. Kuro SJ
par(new = TRUE)
plot(as.matrix(dist1)[10:16, 20:31],
     as.matrix(fst.dist)[10:16, 20:31],
     col = "red", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 2000), cex = 1.5)

# Kuro SJ vs. Kuro SJ
par(new = TRUE)
plot(as.matrix(dist1)[20:31, 20:31],
     as.matrix(fst.dist)[20:31, 20:31],
     col = "red", pch = 22,
     ylim = c(0.00, 0.30), xlim = c(0, 2000), cex = 1.5)

dev.off()


########################################
# 8. Environmental PCA
########################################

env_data <- read.table("env.txt", header = TRUE)

env_pca <- prcomp(env_data, scale. = TRUE)

# Inspect PCA variance explained.
print(summary(env_pca))

pca_scores <- as.data.frame(env_pca$x)

pop_names <- c(
    "PO_1", "PO_2", "PO_3", "PO_4", "PO_5", "PO_6", "PO_7", "PO_8",
    "PO_9", "PO_10", "PO_11", "PO_12", "PO_13", "PO_14", "PO_15", "PO_16",
    "SJ_1", "SJ_2", "SJ_3", "SJ_4", "SJ_5", "SJ_6", "SJ_7", "SJ_8",
    "SJ_9", "SJ_10", "SJ_11", "SJ_12", "SJ_13", "SJ_14", "SJ_15"
)

pca_scores$pop <- factor(pop_names, levels = pop_names)
rownames(pca_scores) <- pop_names

colors <- colorRampPalette(
    c("light blue", "blue", "purple", "green", "orange", "red")
)

col_list <- c(colors(16), colors(15))


# PC1-PC2 plot.
g <- ggplot(
    pca_scores,
    aes(x = PC1, y = PC2, color = pop)
) +
    geom_point(size = 4) +
    scale_color_manual(values = col_list) +
    ggrepel::geom_label_repel(
        aes(label = rownames(pca_scores))
    )

print(g)


# PC scores.
PC1score <- env_pca$x[, 1]
PC2score <- env_pca$x[, 2]

pc1_data <- data.frame(
    pop = pca_scores$pop,
    PC1score = PC1score
)

pc2_data <- data.frame(
    pop = pca_scores$pop,
    PC2score = PC2score
)

x <- 1:31


# PC1 plot.
pdf("PC1_env.pdf", width = 10, height = 8)

par(mar = c(10, 5, 5, 5))
par(bg = "transparent")

beeswarm(
    PC1score ~ pop,
    data = pc1_data,
    las = 2,
    col = col_list,
    pch = 16,
    cex = 2.5,
    xlab = "",
    ylab = "PC1 score"
)

abline(v = x, col = "grey", lty = 3)

par(new = TRUE)

beeswarm(
    PC1score ~ pop,
    data = pc1_data,
    las = 2,
    col = col_list,
    pch = 16,
    cex = 2.5,
    xlab = "",
    ylab = ""
)

dev.off()


# PC2 plot.
pdf("PC2_env.pdf", width = 13, height = 5)

par(mar = c(10, 5, 5, 5))
par(bg = "transparent")

beeswarm(
    PC2score ~ pop,
    data = pc2_data,
    las = 2,
    col = col_list,
    pch = 16,
    xlab = "",
    ylab = "PC2 score"
)

abline(v = x, col = "grey", lty = 3)

par(new = TRUE)

beeswarm(
    PC2score ~ pop,
    data = pc2_data,
    las = 2,
    col = col_list,
    pch = 16,
    xlab = "",
    ylab = ""
)

dev.off()


########################################
# 9. Mantel test: isolation by environment
########################################

# Environmental distance based on PC1 scores.
all_env_dist <- dist(PC1score, method = "euclidean")

ibe <- vegan::mantel(fst.dist, all_env_dist)

sink("env_Mantel_result.txt")
print(ibe)
sink()


########################################
# 10. IBE plots
########################################

# PO vs. SJ comparison.
pdf("IBE_FST_GRAS.pdf", width = 10, height = 10)

par(mar = c(10, 5, 10, 5))

plot(
    all_env_dist,
    fst.dist,
    ylim = c(0.00, 0.30),
    xlim = c(0, 16)
)

# PO vs. SJ
par(new = TRUE)
plot(
    as.matrix(all_env_dist)[17:31, 1:16],
    as.matrix(fst.dist)[17:31, 1:16],
    col = "blue",
    pch = 16,
    ylim = c(0.00, 0.30),
    xlim = c(0, 16)
)

# PO vs. PO
par(new = TRUE)
plot(
    as.matrix(all_env_dist)[1:16, 1:16],
    as.matrix(fst.dist)[1:16, 1:16],
    col = "orange",
    pch = 16,
    ylim = c(0.00, 0.30),
    xlim = c(0, 16)
)

# SJ vs. SJ
par(new = TRUE)
plot(
    as.matrix(all_env_dist)[17:31, 17:31],
    as.matrix(fst.dist)[17:31, 17:31],
    col = "green",
    pch = 16,
    ylim = c(0.00, 0.30),
    xlim = c(0, 16)
)

dev.off()


# Ezo vs. Kuro comparison.
pdf("IBE_FST_GRAS_R1.pdf", width = 10, height = 10)

par(mar = c(10, 5, 10, 5))

# Ezo vs. Kuro
plot(as.matrix(all_env_dist)[1:9, 10:16],
     as.matrix(fst.dist)[1:9, 10:16],
     col = "grey", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 16), cex = 1.5)

par(new = TRUE)
plot(as.matrix(all_env_dist)[1:9, 20:31],
     as.matrix(fst.dist)[1:9, 20:31],
     col = "grey", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 16), cex = 1.5)

par(new = TRUE)
plot(as.matrix(all_env_dist)[10:16, 17:19],
     as.matrix(fst.dist)[10:16, 17:19],
     col = "grey", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 16), cex = 1.5)

par(new = TRUE)
plot(as.matrix(all_env_dist)[17:19, 20:31],
     as.matrix(fst.dist)[17:19, 20:31],
     col = "grey", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 16), cex = 1.5)

# Ezo
par(new = TRUE)
plot(as.matrix(all_env_dist)[1:9, 1:9],
     as.matrix(fst.dist)[1:9, 1:9],
     col = "blue", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 16), cex = 1.5)

par(new = TRUE)
plot(as.matrix(all_env_dist)[1:9, 17:19],
     as.matrix(fst.dist)[1:9, 17:19],
     col = "blue", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 16), cex = 1.5)

par(new = TRUE)
plot(as.matrix(all_env_dist)[17:19, 17:19],
     as.matrix(fst.dist)[17:19, 17:19],
     col = "blue", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 16), cex = 1.5)

# Kuro
par(new = TRUE)
plot(as.matrix(all_env_dist)[10:16, 10:16],
     as.matrix(fst.dist)[10:16, 10:16],
     col = "red", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 16), cex = 1.5)

par(new = TRUE)
plot(as.matrix(all_env_dist)[10:16, 20:31],
     as.matrix(fst.dist)[10:16, 20:31],
     col = "red", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 16), cex = 1.5)

par(new = TRUE)
plot(as.matrix(all_env_dist)[20:31, 20:31],
     as.matrix(fst.dist)[20:31, 20:31],
     col = "red", pch = 21,
     ylim = c(0.00, 0.30), xlim = c(0, 16), cex = 1.5)

dev.off()


########################################
# 11. Partial Mantel test
########################################

# Association between genetic and environmental distances
# while controlling for geographic distance.
partial <- vegan::mantel.partial(
    fst.dist,
    all_env_dist,
    dist1
)

sink("partial_Mantel_result.txt")
print(partial)
sink()
