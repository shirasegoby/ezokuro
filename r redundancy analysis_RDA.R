

library(vegan)
library(dplyr)



##################################################
gen <- read.csv('aba_geno.csv', header = TRUE, row.names=1)
env <- read.csv('aba_env.csv', header = TRUE, row.names=1)
env <-
  env %>%
  mutate(lg.N.q1 = log10(N.q1)) %>%
  mutate(lg.N.q2 = log10(N.q2)) %>%
  mutate(lg.N.q3 = log10(N.q3)) %>%
  mutate(lg.N.q4 = log10(N.q4)) %>%
  mutate(lg.P.q1 = log10(P.q1)) %>%
  mutate(lg.P.q2 = log10(P.q2)) %>%
  mutate(lg.P.q3 = log10(P.q3)) %>%
  mutate(lg.P.q4 = log10(P.q4)) %>%
  mutate(lg.Si.q1 = log10(Si.q1)) %>%
  mutate(lg.Si.q2 = log10(Si.q2)) %>%
  mutate(lg.Si.q3 = log10(Si.q3)) %>%
  mutate(lg.Si.q4 = log10(Si.q4))

env <- env[, c(1:9,30:45)] # starting with a total of 25 variables (fetch, temp q1:4, sali q1:4, pH q1:4, N q1:4, P q1:4, Si q1:4)
dim(gen) # [1] 154 448 (= 68,992 data in the matrix)
sum(is.na(gen)) # 1265 NAs in the matrix (~1.8% missing data)
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) #imputing missing data (NAs) using the most common genotype at each SNP across all individuals.


# rda with all the environmental variables
aba.rda <- rda(gen.imp ~ ., data=env, scale=T)
aba.rda

RsquareAdj(aba.rda)
summary(eigenvals(aba.rda, model = "constrained"))
screeplot(aba.rda)


# step.forward selection
rda.aba1 <- rda(gen.imp~., data = env, scale=T)
rda.aba0 <- rda(gen.imp ~ 1, data = env, scale=T)

step.forwardL <- 
  ordistep(rda.aba0, 
           scope = formula(rda.aba1), 
           direction = "forward", 
           permutations = how(nperm = 999)
  )

RsquareAdj(step.forwardL)
step.forwardL$anova
plot(step.forwardL)
anova(step.forwardL)
coef(step.forwardL)
summary(step.forwardL)


# rda with best env variables
env.best <- env[, c(1,5,13,16,20,25)] # based on the forward selection
colnames(env.best) <- c("FETCH","TEMP.Q4", "pH.Q4", "N.Q3", "P.Q3", "Si.Q4") 
vifL <- diag(solve(cor(env.best))) # checking VIF values for the six variables
vifL

rda.best <- rda(gen.imp~., data = env.best, scale=T)
anova(rda.best, permutations = how(nperm = 999))
RsquareAdj(rda.best)
plot(rda.best)
anova(rda.best)
summary(rda.best)


# Plot the rda across groups
grp <- read.csv('aba_group.csv', header = TRUE, row.names=1)
grp$Group <- as.factor(grp$Group)
levels(grp$Group) <- c("A","B","C","D")
aba <- grp$Group
bg <- c("#1f78b4","#a6cee3","#ff7f00","#ffff33")

op<-par(mfrow=c(1,1), mar=c(4.5, 4.5, 2, 1))
plot(rda.best, scaling = 3, type = "n", col = "black")
points(rda.best, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[aba])
legend("topleft", legend=levels(aba), bty="n", col="gray32", pch=21, cex=0.9, pt.bg=bg)
mtext("(95.3 % of fitted, 53.6 % of total variation)", side = 1, adj = 0.5, line = 1.9, cex = 0.8)
mtext("(2.3 % of fitted, 1.3 % of total variation)", side = 2, adj = 0.5, line = 2.1, cex = 0.8)
ef <- envfit (rda.best, env.best, scaling = 3, permutations = 0)
plot (ef, cex = 0.8)
par(op)
dev.off()

