 # Load (and install, if needed) packages and convenience functions
source("scripts/setup.R")
source("scripts/headers.R")

####################################
# Load Clean Data ##################
####################################

alpha.div <- readRDS("clean_data/alphaDiversity.rds")
beta.div <- readRDS("clean_data/betaDiversity.rds")
env <- readRDS("raw_data/climate_data.RDS")
tmp <- readRDS("raw_data/cru/cru-tmp.RDS") 

alpha.div$Temp <- -1*alpha.div$PC1
alpha.div <- as.data.frame(scale(alpha.div))
beta.div$Temp <- -1*beta.div$pc1
beta.div <- as.data.frame(scale(beta.div))

###################################
#### Correlations among metrics ###
###################################

## Table 1
# note: does not indicate p-value

ab.corr <- cbind(as.data.frame(c("Alpha Diversity", "", "",  "Beta Diversity", "", "")),
                 as.data.frame(c("tax-phy", "tax-fun", "phy-fun", "tax-PCDp", "tax-fun", "PCDp-fun")),
                 as.data.frame(c(cor.test(alpha.div$plant.tax, alpha.div$plant.phy)$estimate, 
                                 cor.test(alpha.div$plant.tax, alpha.div$plant.fun)$estimate,
                                 cor.test(alpha.div$plant.phy, alpha.div$plant.fun)$estimate,
                                 cor.test(beta.div$plant.tax, beta.div$plant.pcdp)$estimate,
                                 cor.test(beta.div$plant.tax, beta.div$plant.fun)$estimate,
                                 cor.test(beta.div$plant.pcdp, beta.div$plant.fun)$estimate)),
                 as.data.frame(c(cor.test(alpha.div$mammal.tax,alpha.div$mammal.phy)$estimate, 
                                cor.test(alpha.div$mammal.tax, alpha.div$mammal.fun)$estimate,
                                cor.test(alpha.div$mammal.phy, alpha.div$mammal.fun)$estimate,
                                cor.test(beta.div$mammal.tax, beta.div$mammal.pcdp)$estimate,
                                cor.test(beta.div$mammal.tax, beta.div$mammal.fun)$estimate,
                                cor.test(beta.div$mammal.pcdp, beta.div$mammal.fun)$estimate)),
                 as.data.frame(c(cor.test(alpha.div$bird.tax,alpha.div$bird.phy)$estimate, 
                                cor.test(alpha.div$bird.tax, alpha.div$bird.fun)$estimate,
                                cor.test(alpha.div$bird.phy, alpha.div$bird.fun)$estimate,
                                cor.test(beta.div$bird.tax, beta.div$bird.pcdp)$estimate,
                                cor.test(beta.div$bird.tax, beta.div$bird.fun)$estimate,
                                cor.test(beta.div$bird.pcdp, beta.div$bird.fun)$estimate)))
names(ab.corr) <- c("", "Comparison", "Plants", "Mammals", "Birds")
print(xtable(ab.corr), include.rownames=FALSE, digits=c(0,0,2,2,2))

# MANTEL TESTS for beta diversity

plant.tax <- convert.dist(beta.div$plant.tax)
plant.pcdp <- convert.dist(beta.div$plant.pcdp)
plant.fun <- convert.dist(beta.div$plant.fun)

mammal.tax <- convert.dist(beta.div$mammal.tax)
mammal.pcdp <- convert.dist(beta.div$mammal.pcdp)
mammal.fun <- convert.dist(beta.div$mammal.fun)

bird.tax <- convert.dist(beta.div$bird.tax)
bird.pcdp <- convert.dist(beta.div$bird.pcdp)
bird.fun <- convert.dist(beta.div$bird.fun)

mantel(plant.tax, plant.pcdp, method= "pearson", permutations = 999)
mantel(plant.tax, plant.fun, method= "pearson", permutations = 999)
mantel(plant.pcdp, plant.fun, method= "pearson", permutations = 999)

mantel(mammal.tax, mammal.pcdp, method= "pearson", permutations = 999, na.rm = TRUE)
mantel(mammal.tax, mammal.fun, method= "pearson", permutations = 999)
mantel(mammal.pcdp, mammal.fun, method= "pearson", permutations = 999, na.rm = TRUE)

mantel(bird.tax, bird.pcdp, method= "pearson", permutations = 999)
mantel(bird.tax, bird.fun, method= "pearson", permutations = 999)
mantel(bird.pcdp, bird.fun, method= "pearson", permutations = 999)

#############################
####  Linear Models of ######
####  Alpha Diversity  ######
#############################

#### all explanatory variables

# Plants
#taxonomic
plant.tax.all <- summary(lm(data=alpha.div, plant.tax ~ mammal.tax + bird.tax + 
                              Temp + PC2 + Elevation))
p.tax <- coefs(plant.tax.all)
#phylogenetic
plant.phy.all <- summary(lm(data=alpha.div, plant.phy ~ mammal.phy + bird.phy + 
                              Temp + PC2 + Elevation))
p.phy <- coefs(plant.phy.all)
#functional
plant.fun.all <- summary(lm(data=alpha.div, plant.fun ~ mammal.fun + bird.fun +  
                              Temp + PC2 + Elevation))
p.fun <-  coefs(plant.fun.all)

# Birds
#taxonomic
bird.tax.all <- summary(lm(data=alpha.div, bird.tax ~ plant.tax + mammal.tax +  
                             Temp + PC2 + Elevation))
b.tax <- coefs(bird.tax.all)
#phylogenetic
bird.phy.all <- summary(lm(data=alpha.div, bird.phy ~ plant.phy + mammal.phy +  
                             Temp + PC2 + Elevation))
b.phy <- coefs(bird.phy.all)
#functional
bird.fun.all <- summary(lm(data=alpha.div, bird.fun ~ plant.fun + mammal.fun +  
                             Temp + PC2 + Elevation))
b.fun <- coefs(bird.fun.all)

# Mammals
#taxonomic
mam.tax.all <- summary(lm(data=alpha.div, mammal.tax ~ plant.tax + bird.tax +  
                            Temp + PC2 + Elevation))
m.tax <- coefs(mam.tax.all)
#phylogenetic
mam.phy.all <- summary(lm(data=alpha.div, mammal.phy ~ plant.phy + bird.phy +  
                            Temp + PC2 + Elevation))
m.phy <- coefs(mam.phy.all)
#functional
mam.fun.all <- summary(lm(data=alpha.div, mammal.fun ~ plant.fun + bird.fun +  
                            Temp + PC2 + Elevation))
m.fun <- coefs(mam.fun.all)


#### only abiotic explanatory variables

plant.tax.ab <- lm(data=alpha.div, plant.tax ~ Temp + PC2 + Elevation)
#phylogenetic
plant.phy.ab <- lm(data=alpha.div, plant.phy ~ Temp + PC2 + Elevation)
#functional
plant.fun.ab <- lm(data=alpha.div, plant.fun ~ Temp + PC2 + Elevation)

# Birds
#taxonomic
bird.tax.ab <- lm(data=alpha.div, bird.tax ~ Temp + PC2 + Elevation)
#phylogenetic
bird.phy.ab <- lm(data=alpha.div, bird.phy ~ Temp + PC2 + Elevation)
#functional
bird.fun.ab <- lm(data=alpha.div, bird.fun ~ Temp + PC2 + Elevation)

# Mammals
#taxonomic
mam.tax.ab <- lm(data=alpha.div, mammal.tax ~ Temp + PC2 + Elevation)
#phylogenetic
mam.phy.ab <- lm(data=alpha.div, mammal.phy ~ Temp + PC2 + Elevation)
#functional
mam.fun.ab <- lm(data=alpha.div, mammal.fun ~ Temp + PC2 + Elevation)

#### only biotic explanatory variables

# Plants
#taxonomic
plant.tax.bi <- summary(lm(data=alpha.div, plant.tax ~ mammal.tax + bird.tax))
#phylogenetic
plant.phy.bi <- summary(lm(data=alpha.div, plant.phy ~ mammal.phy + bird.phy))
#functional
plant.fun.bi <- summary(lm(data=alpha.div, plant.fun ~ mammal.fun + bird.fun))

# Birds
#taxonomic
bird.tax.bi <- summary(lm(data=alpha.div, bird.tax ~ plant.tax + mammal.tax))
#phylogenetic
bird.phy.bi <- summary(lm(data=alpha.div, bird.phy ~ plant.phy + mammal.phy))
#functional
bird.fun.bi <- summary(lm(data=alpha.div, bird.fun ~ plant.fun + mammal.fun))

# Mammals
#taxonomic
mam.tax.bi <- summary(lm(data=alpha.div, mammal.tax ~ plant.tax + bird.tax))
#phylogenetic
mam.phy.bi <- summary(lm(data=alpha.div, mammal.phy ~ plant.phy + bird.phy))
#functional
mam.fun.bi <- summary(lm(data=alpha.div, mammal.fun ~ plant.fun + bird.fun))

#### Abiotic model residuals ~ biotic explanatory variables

# Plants
#taxonomic
plant.tax.resid <- summary(lm(plant.tax.ab$residuals~alpha.div$mammal.tax + alpha.div$bird.tax))
#phylogenetic
plant.phy.resid <- summary(lm(plant.phy.ab$residuals~alpha.div$mammal.phy + alpha.div$bird.phy))
#functional
plant.fun.resid <- summary(lm(plant.fun.ab$residuals~alpha.div$mammal.fun + alpha.div$bird.fun))

# Birds
#taxonomic
bird.tax.resid <- summary(lm(bird.tax.ab$residuals~alpha.div$plant.tax + alpha.div$mammal.tax))
#phylogenetic
bird.phy.resid <- summary(lm(bird.phy.ab$residuals~alpha.div$plant.phy + alpha.div$mammal.phy))
#functional
bird.fun.resid <- summary(lm(bird.fun.ab$residuals~alpha.div$plant.fun + alpha.div$mammal.fun))

# Mammals
#taxonomic
mam.tax.resid <- summary(lm(mam.tax.ab$residuals~alpha.div$plant.tax + alpha.div$bird.tax))
#phylogenetic -- !! 
# this does not work because there are so many NAs in the mammal ses.mntd...
mam.resid<- as.data.frame(mam.phy.ab$residuals); names(mam.resid) <- "residuals"
resid.dat <- alpha.div[row.names(alpha.div) %in% rownames(mam.resid), ]
mam.phy.resid <- summary(lm(mam.resid$residuals~resid.dat$plant.phy + resid.dat$bird.phy))
#functional
mam.fun.resid <- summary(lm(mam.fun.ab$residuals~alpha.div$plant.fun + alpha.div$bird.fun))


## Table 2: comparison of R2 values from alpha diversity models (full, abiotic only, biotic only)
alpha.r2 <- cbind(as.data.frame(c("Taxonomic", "Taxonomic", "Taxonomic", "Phylogenetic", "Phylogenetic",
                                  "Phylogenetic", "Functional", "Functional", "Functional")),
                  as.data.frame(c("Full Model", "Abiotic", "Biotic", 
                                  "Full Model", "Abiotic", "Biotic",
                                  "Full Model", "Abiotic", "Biotic")),
                  as.data.frame(c(plant.tax.all$r.squared, summary(plant.tax.ab)$r.squared, plant.tax.bi$r.squared, 
                                  plant.phy.all$r.squared, summary(plant.phy.ab)$r.squared,
                                  plant.phy.bi$r.squared, plant.fun.all$r.squared,
                                  summary(plant.fun.ab)$r.squared, plant.fun.bi$r.squared)),
                  as.data.frame(c(bird.tax.all$r.squared, summary(bird.tax.ab)$r.squared, bird.tax.bi$r.squared, 
                                  bird.phy.all$r.squared, summary(bird.phy.ab)$r.squared,
                                  bird.phy.bi$r.squared, bird.fun.all$r.squared,
                                  summary(bird.fun.ab)$r.squared, bird.fun.bi$r.squared)),
                  as.data.frame(c(mam.tax.all$r.squared, summary(mam.tax.ab)$r.squared, mam.tax.bi$r.squared, 
                                  mam.phy.all$r.squared, summary(mam.phy.ab)$r.squared,
                                  mam.phy.bi$r.squared, mam.fun.all$r.squared,
                                  summary(mam.fun.ab)$r.squared, mam.fun.bi$r.squared)))
names(alpha.r2) <- c("", "", "Plants", "Birds", "Mammals")
print(xtable(alpha.r2), include.rownames=FALSE)


# make a bar figure
names(alpha.r2) <- c("Type", "Model", "Plants", "Birds", "Mammals")
alpha.r2$Model <- factor(alpha.r2$Model, levels = c("Full Model", "Abiotic", "Biotic"))
alpha.r2$Type <- factor(alpha.r2$Type, levels = c("Taxonomic", "Phylogenetic", "Functional"))
alpha.r2.plant <- alpha.r2[,c(1:3)]
alpha.r2.mammal <- alpha.r2[,c(1:2,5)]
alpha.r2.bird <- alpha.r2[,c(1:2,4)]

cbPalette <- c("#999999", "#56B4E9", "#FF6666", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf("figures/save_out/modelsummaries_plants.pdf", width=6, height=4)
ggplot(alpha.r2.plant, aes(fill=Model, y=Plants, x=Type)) + ylim(0,1) +
  geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=cbPalette) +
  theme(panel.background = element_blank())
dev.off()

pdf("figures/save_out/modelsummaries_mammals.pdf", width=6, height=4)
ggplot(alpha.r2.mammal, aes(fill=Model, y=Mammals, x=Type)) + ylim(0,1) +
  geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=cbPalette) + 
  theme(panel.background = element_blank())
dev.off()

pdf("figures/save_out/modelsummaries_birds.pdf", width=6, height=4)
ggplot(alpha.r2.bird, aes(fill=Model, y=Birds, x=Type)) + ylim(0,1) +
  geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=cbPalette) +
  theme(panel.background = element_blank())
dev.off()

# Same table but with the residual models (abiotic only ~ biotic only models) for the Supplement
# that mammal phylogenetic residual model keeps this from working for now...
alpha.r2.supp <- cbind(as.data.frame(c("Taxonomic", "", "", "", "Phylogenetic", "", "", "", "Functional", "", "", "")),
                  as.data.frame(c("Full Model", "Abiotic", "Biotic",  "Abiotic~Biotic", 
                                  "Full Model", "Abiotic", "Biotic",  "Abiotic~Biotic",
                                  "Full Model", "Abiotic", "Biotic",  "Abiotic~Biotic")),
                  as.data.frame(c(plant.tax.all$r.squared, summary(plant.tax.ab)$r.squared, plant.tax.bi$r.squared, 
                                  plant.tax.resid$r.squared, plant.phy.all$r.squared, summary(plant.phy.ab)$r.squared,
                                  plant.phy.bi$r.squared, plant.phy.resid$r.squared, plant.fun.all$r.squared,
                                  summary(plant.fun.ab)$r.squared, plant.fun.bi$r.squared, plant.fun.resid$r.squared)),
                  as.data.frame(c(bird.tax.all$r.squared, summary(bird.tax.ab)$r.squared, bird.tax.bi$r.squared, 
                                  bird.tax.resid$r.squared, bird.phy.all$r.squared, summary(bird.phy.ab)$r.squared,
                                  bird.phy.bi$r.squared, bird.phy.resid$r.squared, bird.fun.all$r.squared,
                                  summary(bird.fun.ab)$r.squared, bird.fun.bi$r.squared, bird.fun.resid$r.squared)),
                  as.data.frame(c(mam.tax.all$r.squared, summary(mam.tax.ab)$r.squared, mam.tax.bi$r.squared, 
                                  mam.tax.resid$r.squared, mam.phy.all$r.squared, summary(mam.phy.ab)$r.squared,
                                  mam.phy.bi$r.squared, mam.phy.resid$r.squared, mam.fun.all$r.squared,
                                  summary(mam.fun.ab)$r.squared, mam.fun.bi$r.squared, mam.fun.resid$r.squared)))
names(alpha.r2.supp) <- c("", "", "Plants", "Birds", "Mammals")
print(xtable(alpha.r2.supp), include.rownames=FALSE)

## Figure 2A, magnitude of individual coefficients

#combine and label all of the alpha coefs
sub <- c(0,0,0)
plant.coef <- as.data.frame(cbind(p.tax$est, p.phy$est, p.fun$est))
plant.coef <- InsertRow(plant.coef, NewRow=sub, RowNum=1)
mam.coef <- as.data.frame(cbind(m.tax$est, m.phy$est, m.fun$est))
mam.coef <- InsertRow(mam.coef, NewRow=sub, RowNum=2)
bird.coef <- as.data.frame(cbind(b.tax$est, b.phy$est, b.fun$est)) 
bird.coef <- InsertRow(bird.coef, NewRow=sub, RowNum=3)
alpha.coef <- cbind(plant.coef, mam.coef, bird.coef)
alpha.coef[is.na(alpha.coef)] = 0
rownames(alpha.coef) <- c("Plant", "Mammal", "Bird", "Temperature", "Precipitation", "Elevation")
colnames(alpha.coef) <- c("Plant Tax", "Plant Phy", "Plant Fun","Mammal Tax", "Mammal Phy", 
                          "Mammal Fun","Bird Tax", "Bird Phy", "Bird Fun")

coef <- as.matrix(alpha.coef)
coef <- t(coef)

# p-values for coefs
sub <-  c(1,1,1)
plant.p <- as.data.frame(cbind(p.tax$p, p.phy$p, p.fun$p))
plant.p <- InsertRow(plant.p, NewRow=sub, RowNum=1)
mam.p <- as.data.frame(cbind(m.tax$p, m.phy$p, m.fun$p))
mam.p <- InsertRow(mam.p, NewRow=sub, RowNum=2)
bird.p <- as.data.frame(cbind(b.tax$p, b.phy$p, b.fun$p)) 
bird.p <- InsertRow(bird.p, NewRow=sub, RowNum=3)
alpha.p <- cbind(plant.p, mam.p, bird.p)
alpha.p[is.na(alpha.p)] = 0
rownames(alpha.p) <- c("Plant", "Mammal", "Bird", "Temperature", "Precipitation", "Elevation")
colnames(alpha.p) <- c("Plant Tax", "Plant Phy", "Plant Fun","Mammal Tax", "Mammal Phy", 
                        "Mammal Fun","Bird Tax", "Bird Phy", "Bird Fun")

p.val <- as.matrix(alpha.p)
p.val <- t(p.val)

png(filename="figures/save_out/alpha_corrplot_allPval.png", 
    units="in", 
    width=7, 
    height=5, 
    res=150)
par(mfrow=c(1,1), mai=c(1,1,1,1), oma=c(1,1,1,1))
col <- colorRampPalette(c("#990000", "#FF9999", "#FFFFFF", "#6699FF", "#000066"))
corrplot(coef, method="circle", is.corr=FALSE, tl.col="black", tl.cex=0.8,
         cl.ratio = 0.4, cl.align = "r", col=col(200), cl.lim = c(-1.5, 1.5), p.mat=p.val,
         sig.level= 0.05, insig="label_sig", pch=".", pch.cex=2)
dev.off()


###### Figure 1. Map of bird species richness (alpha diversity) 

# linear models of bird richness with Temp and plant richness
par(mfrow=c(1,1))

pdf("figures/save_out/birdPhy_Temp.pdf", width=6, height=4)
plot(alpha.div$Temp, alpha.div$bird.phy, axes=FALSE, col=rgb(1,1,1,0.5), xlim=c(-2.5,3), ylim=c(-3,3),
     xlab="Temperature", ylab="Bird Phylogentic Diversity")
axis(1)
axis(2)
text(bird.phy~Temp, labels=rownames(alpha.div),data=alpha.div, cex=0.7, font=1)
abline(lm(alpha.div$bird.phy~alpha.div$Temp), col="red")
dev.off()

pdf("figures/save_out/birdPhy_plantPhy.pdf", width=6, height=4)
plot(alpha.div$plant.phy, alpha.div$bird.phy, axes=FALSE, col=rgb(1,1,1,0.5), xlim=c(-2,3), ylim=c(-3,3),
     xlab="Plant Phylogenetic Diversity", ylab="Bird Phylogenetic Diversity")
axis(1)
axis(2)
text(bird.phy~plant.phy, labels=rownames(alpha.div),data=alpha.div, cex=0.7, font=1)
abline(lm(alpha.div$bird.phy~alpha.div$plant.phy), col="red")
dev.off()

###### the map part of the map

# get USA lat/lon, separate contiguous
usa.dat <- getData('GADM', country='USA', level=1)
pr.dat <- getData('GADM', country='Puerto Rico', level=1)
simple.usa <- gSimplify(usa.dat, tol=0.1)
simple.pr <- gSimplify(pr.dat, tol=0.1, topologyPreserve=TRUE)
simple.usa <- SpatialPolygonsDataFrame(simple.usa, data=usa.dat@data)
simple.pr <- SpatialPolygonsDataFrame(simple.pr, data=pr.dat@data)
plot(simple.usa)
plot(simple.pr)
usa_contig <- simple.usa[-which(simple.usa$NAME_1 %in% c("Alaska", "Hawaii")),]
pr <- simple.pr

# estimate average 30 year temp across map
stackavg <- calc(x=tmp[[89:118]], fun=mean)
temp_avg_contig <- mask(stackavg, usa_contig)

temp_avg_pr <- mask(stackavg, pr)

# make the figure
pdf("figures/contig_NEONtemps.pdf", width=6, height=4)
cuts=c(-5, 0, 5, 10, 15, 20, 25, 30)
pal <- colorRampPalette(c("skyblue", "white", "red"))
plot(temp_avg_contig, ylim=c(20,60), xlim=c(-140,-60), 
     breaks=cuts, col = pal(length(cuts)), bty="n", box=FALSE)

# add the states
plot(usa_contig, add=TRUE, border="gray")

# add the neon sites
neon.coords <- env[,c(5,4)]
neon.coords <- SpatialPoints(neon.coords)
points(neon.coords, cex=0.5, pch=20)
text(neon.coords, labels=env$site, cex=0.3)
dev.off()

# now just alaska
AK <- simple.usa[which(simple.usa$NAME_1 == "Alaska"),]
temp_ak <- mask(stackavg, AK)

pdf("figures/AK_NEONtemps.pdf", width=6, height=4)
cuts=c(-15,-10,-5,0,5,10,15,20,25,30)
pal <- colorRampPalette(c("skyblue", "white", "red"))
plot(temp_ak, xlim=c(-200,-100), ylim=c(40, 80), 
     breaks=cuts, col = pal(length(cuts)), bty="n", box=FALSE)
plot(AK, add=TRUE, border="gray")
points(neon.coords, cex-0.5, pch=20)
text(neon.coords, labels=env$site, cex=-0.3)
dev.off()

# note: still missing one Puerto Rico site on this map! 
# make the PR figure
pdf("figures/PR_NEONtemps.pdf", width=6, height=4)
cuts=c(-5, 0, 5, 10, 15, 20, 25, 30)
pal <- colorRampPalette(c("skyblue", "white", "red"))
plot(temp_avg_pr, ylim=c(17, 19), xlim=c(-68,-64), 
     breaks=cuts, col = pal(length(cuts)), bty="n", box=FALSE)

# add the states
plot(pr, add=TRUE, border="gray")

# add the neon sites
neon.coords <- env[,c(5,4)]
neon.coords <- SpatialPoints(neon.coords)
points(neon.coords, cex=0.5, pch=20)
text(neon.coords, labels=env$site, cex=0.6)
dev.off()

#####################################
####  Quantile Regression Models ####
####  of Beta diversity      ########
#####################################

#### all explanatory variables

# Plants
#taxonomic
plant.tax.all <- summary.rq(rq(data=beta.div, plant.tax ~ mammal.tax + bird.tax + 
                                 temp + pc2 + elev + dist), se="nid")
bp.tax <- b.coefs(plant.tax.all)
#phylogenetic
plant.phy.all <- summary(rq(data=beta.div, plant.pcdp ~ mammal.pcdp + bird.pcdp + 
                              temp + pc2 + elev + dist), se="nid")
bp.phy <- b.coefs(plant.phy.all)
#functional
plant.fun.all <- summary(rq(data=beta.div, plant.fun ~ mammal.fun + bird.fun + 
                              temp + pc2 + elev + dist), se="nid")
bp.fun <- b.coefs(plant.fun.all)

# Birds
#taxonomic
bird.tax.all <- summary(rq(data=beta.div, bird.tax ~ plant.tax + mammal.tax +
                             temp + pc2 + elev + dist), se="nid")
bb.tax <- b.coefs(bird.tax.all)
#phylogenetic
bird.phy.all <- summary(rq(data=beta.div, bird.pcdp ~ plant.pcdp + mammal.pcdp + 
                             temp + pc2 + elev + dist), se="nid")

bb.phy <- b.coefs(bird.phy.all)
#functional
bird.fun.all <- summary(rq(data=beta.div, bird.fun ~ plant.fun + mammal.fun + 
                             temp + pc2 + elev + dist), se="nid")
bb.fun <- b.coefs(bird.fun.all)

# Mammals
#taxonomic
mam.tax.all <- summary(rq(data=beta.div, mammal.tax ~ plant.tax + bird.tax + 
                            temp + pc2 + elev + dist), se="nid")
bm.tax <- b.coefs(mam.tax.all)

#phylogenetic
mam.phy.all <- summary(rq(data=beta.div, mammal.pcdp ~ plant.pcdp + bird.pcdp + 
                            temp + pc2 + elev + dist), se="nid")
bm.phy <- b.coefs(mam.phy.all)
#functional
mam.fun.all <- summary(rq(data=beta.div, mammal.fun ~ plant.fun + bird.fun + 
                            temp + pc2 + elev + dist), se="nid")
bm.fun <- b.coefs(mam.fun.all)

#### only abiotic explanatory variables

plant.tax.ab <- rq(data=beta.div, plant.tax ~ temp + pc2 + elev + dist)
#phylogenetic
plant.phy.ab <- rq(data=beta.div, plant.pcdp ~ temp + pc2 + elev + dist)
#functional
plant.fun.ab <- rq(data=beta.div, plant.fun ~ temp + pc2 + elev + dist)

# Birds
#taxonomic
bird.tax.ab <- rq(data=beta.div,bird.tax ~ temp + pc2 + elev + dist)
#phylogenetic
bird.phy.ab <- rq(data=beta.div, bird.pcdp ~ temp + pc2 + elev + dist)
#functional
bird.fun.ab <- rq(data=beta.div, bird.fun ~ temp + pc2 + elev + dist)

# Mammals
#taxonomic
mam.tax.ab <- rq(data=beta.div, mammal.tax ~ temp + pc2 + elev + dist)
#phylogenetic
mam.phy.ab <- rq(data=beta.div, mammal.pcdp ~ temp + pc2 + elev + dist)
#functional
mam.fun.ab <- rq(data=beta.div, mammal.fun ~ temp + pc2 + elev + dist)

#### only biotic explanatory variables

# Plants
#taxonomic
plant.tax.bi <- summary(rq(data=beta.div, plant.tax ~ mammal.tax + bird.tax), se="nid")
#phylogenetic
plant.phy.bi <- summary(rq(data=beta.div, plant.pcdp ~ mammal.pcdp + bird.pcdp), se="nid")
#functional
plant.fun.bi <- summary(rq(data=beta.div, plant.fun ~ mammal.fun + bird.fun), se="nid")

# Birds
#taxonomic
bird.tax.bi <- summary(rq(data=beta.div, bird.tax ~ plant.tax + mammal.tax), se="nid")
#phylogenetic
bird.phy.bi <- summary(rq(data=beta.div, bird.pcdp ~ plant.pcdp + mammal.pcdp), se="nid")
#functional
bird.fun.bi <- summary(rq(data=beta.div, bird.fun ~ plant.fun + mammal.fun), se="nid")

# Mammals
#taxonomic
mam.tax.bi <- summary(rq(data=beta.div, mammal.tax ~ plant.tax + bird.tax), se="nid")
#phylogenetic
mam.phy.bi <- summary(rq(data=beta.div, mammal.pcdp ~ plant.pcdp + bird.pcdp), se="nid")
#functional
mam.fun.bi <- summary(rq(data=beta.div, mammal.fun ~ plant.fun + bird.fun), se="nid")

#### Residuals of abiotic models ~ biotic explanatory variables

# Plants
#taxonomic
plant.tax.resid <- summary(rq(plant.tax.ab$residuals~beta.div$mammal.tax + beta.div$bird.tax), se="nid")
#phylogenetic
plant.phy.resid <- summary(rq(plant.phy.ab$residuals~beta.div$mammal.pcdp + beta.div$bird.pcdp), se="nid")
#functional
plant.fun.resid <- summary(rq(plant.fun.ab$residuals~beta.div$mammal.fun + beta.div$bird.fun), se="nid")

# Birds
#taxonomic
bird.tax.resid <- summary(rq(bird.tax.ab$residuals~beta.div$plant.tax + beta.div$mammal.tax), se="nid")
#phylogenetic
bird.phy.resid <- summary(rq(bird.phy.ab$residuals~beta.div$plant.pcdp + beta.div$mammal.pcdp), se="nid")
#functional
bird.fun.resid <- summary(rq(bird.fun.ab$residuals~beta.div$plant.fun + beta.div$mammal.fun), se="nid")

# Mammals
#taxonomic
mam.tax.resid <- summary(rq(mam.tax.ab$residuals~beta.div$plant.tax + beta.div$bird.tax), se="nid")
#phylogenetic
mam.resid<- as.data.frame(mam.phy.ab$residuals); names(mam.resid) <- "residuals"
resid.dat <- beta.div[row.names(beta.div) %in% rownames(mam.resid), ]
mam.phy.resid <- summary(rq(mam.resid$residuals~resid.dat$plant.pcdp + resid.dat$bird.pcdp), se="nid")
#functional
mam.fun.resid <- summary(rq(mam.fun.ab$residuals~beta.div$plant.fun + beta.div$bird.fun), se="nid")

#### Figure 2B, magnitude of individual coefficients

#combine and label all beta coefs from full models
sub <- c(0,0,0)
bplant.coef <- as.data.frame(cbind(bp.tax$est, bp.phy$est, bp.fun$est))
bplant.coef <- InsertRow(bplant.coef, NewRow=sub, RowNum=1)
bmam.coef <- as.data.frame(cbind(bm.tax$est, bm.phy$est, bm.fun$est))
bmam.coef <- InsertRow(bmam.coef, NewRow=sub, RowNum=2)
bbird.coef <- as.data.frame(cbind(bb.tax$est, bb.phy$est, bb.fun$est)) 
bbird.coef <- InsertRow(bbird.coef, NewRow=sub, RowNum=3)

beta.coef <- cbind(bplant.coef, bmam.coef, bbird.coef)
beta.coef[is.na(beta.coef)] = 0
rownames(beta.coef) <- c("Plant", "Mammal", "Bird", "Temperature", "Precipitation", "Elevation", "Distance")
colnames(beta.coef) <- c("Plant Tax", "Plant Phy", "Plant Fun","Mammal Tax", "Mammal Phy", 
                          "Mammal Fun","Bird Tax", "Bird Phy", "Bird Fun")

beta.coef <- as.matrix(beta.coef)
beta.coef <- t(beta.coef)

# now the p-values
sub <- c(1,1,1)
bplant.p <- as.data.frame(cbind(bp.tax$p, bp.phy$p, bp.fun$p))
bplant.p <- InsertRow(bplant.p, NewRow=sub, RowNum=1)
bmam.p <- as.data.frame(cbind(bm.tax$p, bm.phy$p, bm.fun$p))
bmam.p <- InsertRow(bmam.p, NewRow=sub, RowNum=2)
bbird.p <- as.data.frame(cbind(bb.tax$p, bb.phy$p, bb.fun$p)) 
bbird.p <- InsertRow(bbird.p, NewRow=sub, RowNum=3)

beta.p <- cbind(bplant.p, bmam.p, bbird.p)
beta.p[is.na(beta.p)] = 1
rownames(beta.p) <- c("Plant", "Mammal", "Bird", "Temperature", "Precipitation", "Elevation", "Distance")
colnames(beta.p) <- c("Plant Tax", "Plant Phy", "Plant Fun","Mammal Tax", "Mammal Phy", 
                         "Mammal Fun","Bird Tax", "Bird Phy", "Bird Fun")

beta.p <- as.matrix(beta.p)
beta.p <- t(beta.p)

#figure

png(filename="figures/save_out/beta_corrplot_allPval_pcdp.png", 
    units="in", 
    width=7, 
    height=5, 
    res=150)
par(mfrow=c(1,1), mai=c(1,1,1,1), oma=c(1,1,1,1))
col <- colorRampPalette(c("#990000", "#FF9999", "#FFFFFF", "#6699FF", "#000066"))
corrplot(beta.coef, method="circle", is.corr=FALSE, tl.col="black", tl.cex=0.8,
         cl.ratio = 0.4, cl.align = "r",col=col(200), cl.lim = c(-1.5, 1.5),  p.mat=beta.p,
         sig.level= 0.05, insig="label_sig", pch=".", pch.cex=2)
dev.off()

#####  Figure 3A and B: 3D scatterplots of beta diversity

# 3A: bird taxonomic diversity ~ elevation and plant taxonomic diversity

padding <- list(
  layout.heights = list(
    top.padding = 0,
    main.key.padding = 0,
    key.axis.padding = 0,
    axis.xlab.padding = 0,
    xlab.key.padding = 0,
    key.sub.padding = 0,
    bottom.padding = 0
  ),
  layout.widths = list(
    left.padding = 1,
    key.ylab.padding = 0,
    ylab.axis.padding = 1,
    axis.key.padding = 0,
    right.padding = 0
  ),
  axis.line=list(col="transparent"),
  clip =list(panel="off")
)

df <- data.frame(x= beta.div[,1],
                 y= beta.div[,15],
                 z= beta.div[,5])
shadow_x <- df
shadow_y <- df
shadow_z <- df
shadow_x$x <- rep(2.2,nrow(df))
shadow_y$y <- rep(4,nrow(df))
shadow_z$z <- rep(-3,nrow(df))
df_shadows <- rbind(df, shadow_x, shadow_y, shadow_z)

png(filename="figures/save_out/BirdTurnover_3Dscatter.png", 
    units="in", 
    width=6, 
    height=6, 
    res=150)
cloud(df_shadows$z ~ df_shadows$x*df_shadows$y,
      pch=c(rep(1, nrow(shadow_x)), rep(19, nrow(shadow_x)*3)),
      col=c(rep("blue", nrow(shadow_x)), rep("gray", nrow(shadow_x)*3)),
      xlab=list(bquote(atop("Plant Spp.","Turnover")), cex=0.7), 
      ylab=list("Elevation", cex=0.7), 
      zlab=list(bquote(atop("Bird Spp.", "Turnover")), cex=0.7),
      par.settings = padding,
      xlim=c(-6.2,2.2), ylim=c(-1.2,4), zlim=c(-3.2,1.5))
dev.off()

# 3B: mammal taxonomic diversity ~ Temp and plant phylosorenson

df <- data.frame(x= beta.div[,1],
                 y= beta.div[,14],
                 z= beta.div[,9])
shadow_x <- df
shadow_y <- df
shadow_z <- df
shadow_x$x <- rep(2.2,nrow(df))
shadow_y$y <- rep(4,nrow(df))
shadow_z$z <- rep(-3,nrow(df))
df_shadows <- rbind(df, shadow_x, shadow_y, shadow_z)

png(filename="figures/save_out/MammalTurnover_3Dscatter.png", 
    units="in", 
    width=6, 
    height=6, 
    res=150)
cloud(df_shadows$z ~ df_shadows$x*df_shadows$y,
      pch=c(rep(1, nrow(shadow_x)), rep(19, nrow(shadow_x)*3)),
      col=c(rep("blue", nrow(shadow_x)), rep("gray", nrow(shadow_x)*3)),
      xlab=list(bquote(atop("Plant Spp.", "Turnover")), cex=0.7), 
                       ylab=list("Precipitation", cex=0.7), 
      zlab=list(bquote(atop("Mammal Spp.", "Turnover")), cex=0.7),
      par.settings = padding,
      xlim=c(-6,2.2), ylim=c(-2,4), zlim=c(-3,1.2))
dev.off()


##### Below is another version with taxonomic beta diversity 
# (made for ESA presentation) 

# Birds

padding <- list(
  layout.heights = list(
    top.padding = 0,
    main.key.padding = 0,
    key.axis.padding = 0,
    axis.xlab.padding = 0,
    xlab.key.padding = 0,
    key.sub.padding = 0,
    bottom.padding = 0
  ),
  layout.widths = list(
    left.padding = 1,
    key.ylab.padding = 0,
    ylab.axis.padding = 1,
    axis.key.padding = 0,
    right.padding = 0
  ),
  axis.line=list(col="transparent"),
  clip =list(panel="off")
)

df <- data.frame(x= beta.div[,1],
                 y= beta.div[,13],
                 z= beta.div[,5])
shadow_x <- df
shadow_y <- df
shadow_z <- df
shadow_x$x <- rep(2.2,nrow(df))
shadow_y$y <- rep(4,nrow(df))
shadow_z$z <- rep(-3,nrow(df))
df_shadows <- rbind(df, shadow_x, shadow_y, shadow_z)

png(filename="figures/BirdTax_3Dscatter_Jul.png", 
    units="in", 
    width=6, 
    height=6, 
    res=150)
cloud(df_shadows$z ~ df_shadows$x*df_shadows$y,
      pch=c(rep(1, nrow(shadow_x)), rep(19, nrow(shadow_x)*3)),
      col=c(rep("blue", nrow(shadow_x)), rep("gray", nrow(shadow_x)*3)),
      xlab=list("Plant Turnover", cex=1), ylab=list("Temperature", cex=1), 
      zlab=list(bquote(atop("Bird",
                            "Turnover")), cex=1),
      par.settings = padding,
      xlim=c(-3.5,2.2), ylim=c(-1.6,4.1), zlim=c(-3,1.6))
dev.off()

# Mammals

df <- data.frame(x= beta.div[,1],
                 y= beta.div[,13],
                 z= beta.div[,9])
shadow_x <- df
shadow_y <- df
shadow_z <- df
shadow_x$x <- rep(2.2,nrow(df))
shadow_y$y <- rep(4,nrow(df))
shadow_z$z <- rep(-3,nrow(df))
df_shadows <- rbind(df, shadow_x, shadow_y, shadow_z)

png(filename="figures/MammalTax_3Dscatter_Jul.png", 
    units="in", 
    width=6, 
    height=6, 
    res=150)
cloud(df_shadows$z ~ df_shadows$x*df_shadows$y,
      pch=c(rep(1, nrow(shadow_x)), rep(19, nrow(shadow_x)*3)),
      col=c(rep("blue", nrow(shadow_x)), rep("gray", nrow(shadow_x)*3)),
      xlab=list("Plant Turnover", cex=1), ylab=list("Temperature", cex=1), 
      zlab=list(bquote(atop("Mammal",
                            "Turnover")), cex=1),
      par.settings = padding,
      xlim=c(-3.5,2.2), ylim=c(-1.6,4), zlim=c(-3,1.6))
dev.off()

#############################################
### Residual beta-div figures for supplement
##############################################
par(mfrow=c(2,3))

# mammal taxonomic
plot(mam.tax.ab$residuals~beta.div$plant.tax, pch=20, cex=0.5, axes=FALSE, 
     xlab="Plant taxonomic turnover", ylab="Mammal abiotic residuals")
axis(1)
axis(2)
abline(lm(mam.tax.ab$residuals~beta.div$plant.tax), lwd-3, col="red")

# mammal phylogenetic
plot(mam.phy.ab$residuals~beta.div$plant.phy, pch=20, cex=0.5, axes=FALSE, 
     xlab="Plant phylogenetic turnover", ylab="Mammal abiotic residuals")
axis(1)
axis(2)
abline(lm(mam.phy.ab$residuals~beta.div$plant.phy), lwd-3, col="red")

# mammal functional
plot(mam.fun.ab$residuals~beta.div$plant.fun, pch=20, cex=0.5, axes=FALSE, 
     xlab="Plant functional turnover", ylab="Mammal abiotic residuals")
axis(1)
axis(2)
abline(lm(mam.fun.ab$residuals~beta.div$plant.fun), lwd-3, col="red")

# bird taxonomic
plot(bird.tax.ab$residuals~beta.div$plant.tax, pch=20, cex=0.5, axes=FALSE, 
     xlab="Plant taxonomic turnover", ylab="Bird abiotic residuals")
axis(1)
axis(2)
abline(lm(bird.tax.ab$residuals~beta.div$plant.tax), lwd-3, col="red")

# bird phylogenetic
plot(bird.phy.ab$residuals~beta.div$plant.phy, pch=20, cex=0.5, axes=FALSE, 
     xlab="Plant phylogenetic turnover", ylab="Bird abiotic residuals")
axis(1)
axis(2)
abline(lm(bird.phy.ab$residuals~beta.div$plant.phy), lwd-3, col="red")

# bird functional
plot(bird.fun.ab$residuals~beta.div$plant.fun, pch=20, cex=0.5, axes=FALSE, 
     xlab="Plant functional turnover", ylab="Bird abiotic residuals")
axis(1)
axis(2)
abline(lm(bird.fun.ab$residuals~beta.div$plant.fun), lwd-3, col="red")

