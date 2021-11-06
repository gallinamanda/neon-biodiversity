# Load (and install, if needed) packages
source("scripts/setup.R")
source("scripts/headers.R")

########################
#### Load in data  #####
########################

# Load assemblage data
plant <- readRDS("raw_data/NEONplant_raw_data.rds")
bird <- readRDS("raw_data/NEONbird_raw_data.rds")
mammal <- readRDS("raw_data/NEONmammal_raw_data.rds")

# Load traits
plant.trait <- readRDS("raw_data/PlantTraits.rds")
mammal.trait <- readRDS("raw_data/MammalTraits.rds")
bird.trait <- readRDS("raw_data/BirdTraits.rds")

# Load phylogenies
bird.trees <- read.nexus("raw_data/FullEricsonBird.nex")
mammal.trees <- read.nexus("raw_data/FullBDMammal.nex")

# For plant trees, read in and make a list
plant_tree_list <- list.files(path="raw_data/zanne_trees")
plant.trees <- list()
# setwd("~/Dropbox/biotic-env/raw_data/zanne_trees") 

for (i in 1:length(plant_tree_list)){
  one.tree <- read.tree(plant_tree_list[i])
  plant.trees <- append(plant.trees, list(one.tree))
}

# Species names to match community data and phylogenies (birds and mammals)
cleanmammalnames <- read.csv("clean_data/Clean_MammalSppNames.csv")
cleanbirdnames <- read.csv("clean_data/Clean_BirdSppNames.csv")
keep <- readRDS("clean_data/CleanPlantNames.RDS")

# Load in climate data
climate <- readRDS("raw_data/climate_data.RDS")

#############################
#### Clean data   ###########
#############################

#### Make a list of sites common to all three taxa
plant.site <- as.data.frame(unique(substr(plant$data$site.id, 0, 4))); names(plant.site) <- c("site")
mammal.site <- as.data.frame(unique(substr(mammal$data$site.id, 0, 4))); names(mammal.site) <- c("site")
bird.site <- as.data.frame(unique(substr(bird$data$site.id, 0, 4))); names(bird.site) <- c("site")
sites <- merge(plant.site, bird.site, by="site")
sites <- merge(sites, mammal.site, by="site")
# remove BARR because mammal n=1
sites <- as.data.frame(sites[c(1,3:39),]); names(sites) <- c("site")

# limit climate data to represented sites
clim <- merge(climate, sites, by="site")
env <- change.rownames(clim, "site")
env <- env[,c(3:13)]

################
#### Plants ####
################

####### Prep Plant Assemblages
Plant <- plant$data
Plant$site <- str_sub(Plant$site.id,1,4)
Plant$plot <- str_sub(Plant$site.id,6,8)
# limit to species in the shared sites
Plant <- merge(Plant, sites, by="site")
keep <- as.data.frame(keep); names(keep) <- "species"
# merge this refined species list with plant assemlage data
p <- merge(Plant, keep, by="species")
p <- p[,c(2,5,3,1,4)]
p$species <- gsub(" ", "_", p$species)

# change from long to wide and replace NAs with 0s (prep for biodiversity calc loop)
p_wide <- spread(p, species, value)
str(p_wide)
p_wide[is.na(p_wide)] <- 0
str(p_wide)
plant.data <- p_wide

# write this out for use in diversity calculations
saveRDS(plant.data, "clean_data/PlantCommunity.RDS")

####### Prep Plant Phylogeny
plant.names <- as.character(unique(p$species))

# now merge
tree.merge <- function(tree, sp.list){
  new.tree <- congeneric.merge(tree, sp.list, split= "_")
  return(new.tree)
}

plant.merge <- lapply(plant.trees, FUN=tree.merge, sp.list=plant.names)
saveRDS(plant.merge, "clean_data/PreppedPlantTrees.RDS")

####### Prep Plant Traits
plant.names <- as.data.frame(plant.names); names(plant.names) <- "species"
plant.height <- merge(plant.trait, plant.names, by="species")
# we have height for 938/2710 species
plant.height <- change.rownames(plant.height, "species")
saveRDS(plant.height, "clean_data/Clean_PlantMaxHeight.RDS")

#################
#### Mammals ####
#################

####### Prep Mammal Assemblages
Mammal <- mammal$data
Mammal$site <- str_sub(Mammal$site.id,1,4)
Mammal$plot <- str_sub(Mammal$site.id,6,8)
# limit to species in the shared sites
Mammal <- merge(Mammal, sites, by="site")
m.spp <- as.data.frame(unique(Mammal$species)); names(m.spp) <- "species"
# remove all genus-only names (sp. and spp.)
keep <- as.data.frame(m.spp[!grepl("\\.", m.spp$species),]); names(keep) <- "species"
keep <- as.data.frame(word(keep$species, 1,2, sep=" ")); names(keep) <- "species"

m <- merge(Mammal, keep, by="species")
m <- m[,c(2,5,3,1,4)]

# now replace species using the names that match the phylogeny
cleanmammalnames$RawSpecies <- as.character(cleanmammalnames$RawSpecies)
cleanmammalnames$CleanSpecies <- as.character(cleanmammalnames$CleanSpecies)
m$cleanspp <- cleanmammalnames$CleanSpecies[match(m$species, cleanmammalnames$RawSpecies)]
m <- m[,c(1,2,3,6,5)]
names(m)[4] <- "species"
m$species <- gsub(" ", "_", m$species)
m <- aggregate(m$value, by=list(m$site, m$plot, m$site.id, m$species), FUN=sum)
names(m) <- c("site", "plot", "site.id", "species", "value")

# check these species names against phylogeny
m.datnames <- as.data.frame(unique(m$species)); names(m.datnames) <- "species"
m.treenames <- as.data.frame(unique(mammal.trees$tree_2136$tip.label))
names(m.treenames) <- "species"
m.treenames$species <- as.character(m.treenames$species)
setdiff(m.treenames$species, m.datnames$species)
# great!

#### Change from long to wide and replace NAs with 0s
m_wide <- spread(m, species, value)
m_wide[is.na(m_wide)] <- 0
mammal.data <- m_wide

# write this out for use in diversity calculations
saveRDS(mammal.data, "clean_data/MammalCommunity.RDS")

####### Prep Mammal Traits
# limit traits to raw mammal names
str(mammal.trait)
trait.names <- unique(keep[1])
trait.names$species <- gsub(" ", "_", trait.names$species)
mammal.traits <- merge(mammal.trait, trait.names, by="species") 
# missing only one species in traits! Ictidomys_tridecemlineatus

# now give these the clean species names
cleanmammalnames$RawSpecies <- gsub(" ", "_", cleanmammalnames$RawSpecies)
cleanmammalnames$CleanSpecies <- gsub(" ", "_", cleanmammalnames$CleanSpecies)
mammal.traits$cleanspp <- cleanmammalnames$CleanSpecies[match(mammal.traits$species, cleanmammalnames$RawSpecies)]
mammal.traits <- mammal.traits[,c(3,2)]
names(mammal.traits)[1] <- "species"

mammal.traits <- change.rownames(mammal.traits, "species")
saveRDS(mammal.traits, "clean_data/Clean_MammalBodyMass.RDS")

################
#### Birds  ####
################

####### Prep Bird Assemblages
# remove non-overlap sites, and genus-level names
Bird <- bird$data
Bird$site <- str_sub(Bird$site.id,1,4)
Bird$plot <- str_sub(Bird$site.id,6,8)
# limit to species in the shared sites
Bird <- merge(Bird, sites, by="site")
b.spp <- as.data.frame(unique(Bird$species)); names(b.spp) <- "species"
# remove all genus-only names (sp. and spp.)
keep <- as.data.frame(b.spp[!grepl("\\.", b.spp$species),]); names(keep) <- "species"

b <- merge(Bird, keep, by="species")
b <- b[,c(2,5,3,1,4)]

# now replace species using the names that match the phylogeny
b$cleanspp <- cleanbirdnames$CleanSpecies[match(b$species, cleanbirdnames$RawSpecies)]
b <- b[,c(1,2,3,6,5)]
names(b)[4] <- "species"
b$species <- gsub(" ", "_", b$species)
b <- aggregate(b$value, by=list(b$site, b$plot, b$site.id, b$species), FUN=sum)
names(b) <- c("site", "plot", "site.id", "species", "value")

# check these species names against phylogeny
b.datnames <- as.data.frame(unique(b$species)); names(b.datnames) <- "species"
b.treenames <- as.data.frame(unique(bird.trees$tree_2075$tip.label))
names(b.treenames) <- "species"
b.treenames$species <- as.character(b.treenames$species)
setdiff(b.treenames$species, b.datnames$species)
# great!

#### Change from long to wide and replace NAs with 0s
b_wide <- spread(b, species, value)
b_wide[is.na(b_wide)] <- 0
bird.data <- b_wide

# write this out for use in diversity calculations
saveRDS(bird.data, "clean_data/BirdCommunity.RDS")

####### Prep Bird Traits
# limit traits to bird names (using species list from cleaned data set)
str(bird.trait)
bird.traits <- merge(bird.trait, b.datnames, by="species")
# this gives us traits for all 367 species!

bird.traits <- change.rownames(bird.traits, "species")
saveRDS(bird.traits, "clean_data/Clean_BirdBodyMass.RDS")

###############################
##### sampling summaries ######
#### (prep for biodiv. calc) ##
###############################

# Mammals
mam <- as.data.frame(unique(mammal$data$site.id)); names(mam) <- "ID"
mam$date <- str_sub(mam$ID,-5,-1)
mam$siteplot <- str_sub(mam$ID,1,8)
mam$site <- str_sub(mam$ID,1,4)
mam$plot <- str_sub(mam$ID,6,8)
mam.ct <- mam %>%
  group_by(site) %>%
  summarise(dates = n_distinct(date), plots = n_distinct(plot))
mam.ct <- as.data.frame(mam.ct)
names(mam.ct) <- c("site", "mam.dates", "mam.plots")

#check how many dates per plot
m.check <- mam %>%
  group_by(siteplot) %>%
  mutate(visits = n_distinct(date))
m.check <- as.data.frame(m.check)

# Plants
plants <- as.data.frame(unique(plant$data$site.id)); names(plants) <- "ID"
plants$date <- str_sub(plants$ID,-5,-1)
plants$month <- str_sub(plants$ID,-5,-4)
plants$siteplot <- str_sub(plants$ID,1,8)
plants$site <- str_sub(plants$ID,1,4)
plants$plot <- str_sub(plants$ID,6,8)
plant.ct <- plants %>%
  group_by(site) %>%
  summarise(dates = n_distinct(date), plots = n_distinct(plot))
plant.ct <- as.data.frame(plant.ct)
names(plant.ct) <- c("site", "plant.dates", "plant.plots")

#check how many dates per plot
p.check <- plants %>%
  group_by(siteplot) %>%
  mutate(visits = n_distinct(date))
p.check <- as.data.frame(p.check)

# Birds
birds <- as.data.frame(unique(bird$data$site.id)); names(birds) <- "ID"
birds$date <- str_sub(birds$ID,-5,-1)
birds$month <- str_sub(birds$ID,-5,-4)
birds$siteplot <- str_sub(birds$ID,1,8)
birds$site <- str_sub(birds$ID,1,4)
birds$plot <- str_sub(birds$ID,6,8)
bird.ct <- birds %>%
  group_by(site) %>%
  summarise(dates = n_distinct(date), plots = n_distinct(plot))
bird.ct <- as.data.frame(bird.ct)
names(bird.ct) <- c("site", "bird.dates", "bird.plots")

# check how many dates per plot
b.check <- birds %>%
  group_by(siteplot) %>%
  mutate(visits = n_distinct(date))
b.check <- as.data.frame(b.check)

# sampling summary
all.ct <- merge(mam.ct, plant.ct, by="site")
all.ct <- merge(all.ct, bird.ct, by="site")
all.ct
sites <- as.data.frame(all.ct[c(1,3:39),1]); names(sites) <- "site"

# minimum site visit counts (across sites)
# Mammals: 2 (TOOL)
# Plants: 10 (STER)
# Birds: 5 (UKFS)
# minimum number of times each plot was visited: 1

###############################################
## alpha and beta diversity calculations ######
###############################################

################################
######  bird diversity   #######
################################
bird.data <- readRDS("clean_data/BirdCommunity.RDS")
bird.trees <- read.nexus("raw_data/FullEricsonBird.nex")
bird.trees <- bird.trees[1:100]
bird.trait <- readRDS("clean_data/Clean_BirdBodyMass.RDS")
bird.tr <- setNames(bird.trait$log.bodymass, rownames(bird.trait))
str(env) # climate/env data from above

# calculate phylogenetic signal
c(b.phylosig100, b.impute100) %<-% .phylotrait100(b, bird.trees, bird.tr)
b.phylosig <- colMeans(b.phylosig100)

# calculate diversity metrics with sampling plots/dates
c(b.alphatax100, b.alphaphy100, b.alphafun100, b.betatax100, 
  b.betaphy100, b.betapcdp100, b.betafun100) %<-%
  .biodiv100(bird.data, bird.trees, bird.trait, env, sample=TRUE, plot=5)

saveRDS(b.alphafun100, "clean_data/biodiv_outputs/b.alphafun100.RDS")
saveRDS(b.alphaphy100, "clean_data/biodiv_outputs/b.alphaphy100.RDS")
saveRDS(b.alphatax100, "clean_data/biodiv_outputs/b.alphatax100.RDS")
saveRDS(b.betatax100, "clean_data/biodiv_outputs/b.betatax100.RDS")
saveRDS(b.betaphy100, "clean_data/biodiv_outputs/b.betaphy100.RDS")
saveRDS(b.betafun100, "clean_data/biodiv_outputs/b.betafun100.RDS")
saveRDS(b.betapcdp100, "clean_data/biodiv_outputs/b.betapcdp100.RDS")

# or just read them in...
#b.alphafun100 <- readRDS("clean_data/biodiv_outputs/b.alphafun100.RDS")
#b.alphaphy100 <- readRDS("clean_data/biodiv_outputs/b.alphaphy100.RDS")
#b.alphatax100 <- readRDS("clean_data/biodiv_outputs/b.alphatax100.RDS")
#b.betatax100 <- readRDS("clean_data/biodiv_outputs/b.betatax100.RDS")
#b.betaphy100 <- readRDS("clean_data/biodiv_outputs/b.betaphy100.RDS")
#b.betafun100 <- readRDS("clean_data/biodiv_outputs/b.betafun100.RDS")
#b.betapcdp100 <- readRDS("clean_data/biodiv_outputs/b.betapcdp100.RDS")

alpha.bird <- cbind(as.data.frame(colMeans(b.alphatax100[1:100,])), as.data.frame(colMeans(b.alphaphy100[1:100,])),
                    as.data.frame(colMeans(b.alphafun100[1:100,])))
names(alpha.bird) <- c("bird.tax", "bird.phy", "bird.fun")

beta.bird <- cbind(as.data.frame(colMeans(b.betatax100[1:100,])), as.data.frame(colMeans(b.betaphy100[1:100,])),
                   as.data.frame(colMeans(b.betapcdp100[1:100,])), as.data.frame(colMeans(b.betafun100[1:100,])))
names(beta.bird) <- c("bird.tax", "bird.phy", "bird.pcdp", "bird.fun")

#############

# calculate diversity metrics without sampling plots/dates

c(b.alphatax100.all, b.alphaphy100.all, b.alphafun100.all, b.betatax100.all, 
  b.betaphy100.all, b.betapcdp100.all, b.betafun100.all) %<-%
  .biodiv100(bird.data, bird.trees, bird.trait, env, sample=FALSE, plot=5)

saveRDS(b.alphafun100.all, "clean_data/biodiv_outputs/b.alphafun100.all.RDS")
saveRDS(b.alphaphy100.all, "clean_data/biodiv_outputs/b.alphaphy100.all.RDS")
saveRDS(b.alphatax100.all, "clean_data/biodiv_outputs/b.alphatax100.all.RDS")
saveRDS(b.betatax100.all, "clean_data/biodiv_outputs/b.betatax100.all.RDS")
saveRDS(b.betaphy100.all, "clean_data/biodiv_outputs/b.betaphy100.all.RDS")
saveRDS(b.betafun100.all, "clean_data/biodiv_outputs/b.betafun100.all.RDS")
saveRDS(b.betapcdp100.all, "clean_data/biodiv_outputs/b.betapcdp100.all.RDS")

# or just read them in...
#b.alphafun100 <- readRDS("clean_data/biodiv_outputs/b.alphafun100.all.RDS")
#b.alphaphy100 <- readRDS("clean_data/biodiv_outputs/b.alphaphy100.all.RDS")
#b.alphatax100 <- readRDS("clean_data/biodiv_outputs/b.alphatax100.all.RDS")
#b.betatax100 <- readRDS("clean_data/biodiv_outputs/b.betatax100.all.RDS")
#b.betaphy100 <- readRDS("clean_data/biodiv_outputs/b.betaphy100.all.RDS")
#b.betafun100 <- readRDS("clean_data/biodiv_outputs/b.betafun100.all.RDS")
#b.betapcdp100 <- readRDS("clean_data/biodiv_outputs/b.betapcdp100.all.RDS")


alpha.bird.all <- cbind(as.data.frame(colMeans(b.alphatax100.all[1:100,])), as.data.frame(colMeans(b.alphaphy100.all[1:100,])),
                    as.data.frame(colMeans(b.alphafun100.all[1:100,])))
names(alpha.bird.all) <- c("bird.tax.all", "bird.phy.all", "bird.fun.all")

beta.bird.all <- cbind(as.data.frame(colMeans(b.betatax100.all[1:100,])), as.data.frame(colMeans(b.betaphy100.all[1:100,])),
                   as.data.frame(colMeans(b.betapcdp100.all[1:100,])), as.data.frame(colMeans(b.betafun100.all[1:100,])))
names(beta.bird.all) <- c("bird.tax.all", "bird.phy.all", "bird.pcdp.all", "bird.fun.all")


##################################
######  mammal diversity   #######
##################################

mammal.data <- readRDS("clean_data/MammalCommunity.RDS")
mammal.trees <- read.nexus("raw_data/FullBDMammal.nex")
mammal.trees <- mammal.trees[1:100]
mammal.trait <- readRDS("clean_data/Clean_MammalBodyMass.RDS")
str(env) # climate/env data from above
mammal.tree.list <- unique(names(mammal.trees))
mammal.tr <- setNames(mammal.trait$log.bodymass, rownames(mammal.trait))

# estimate phylogenetic signal, impute trait for one species
c(m.phylosig100, m.impute100) %<-% .phylotrait100(m, mammal.trees, mammal.tr)
m.phylosig <- colMeans(m.phylosig100)

# average of imputed trait
m.impute <- colMeans(m.impute100)
m.impute <- as.data.frame(m.impute); names(m.impute) <- "log.bodymass"

# calculate diversity metrics with sampling plots/dates
c(m.alphatax100, m.alphaphy100, m.alphafun100, m.betatax100, 
  m.betaphy100, m.betapcdp100, m.betafun100) %<-%
  .biodiv100(mammal.data, mammal.trees, m.impute, env, sample=TRUE, plot=2)

saveRDS(m.alphafun100, "clean_data/biodiv_outputs/m.alphafun100.RDS")
saveRDS(m.alphaphy100, "clean_data/biodiv_outputs/m.alphaphy100.RDS")
saveRDS(m.alphatax100, "clean_data/biodiv_outputs/m.alphatax100.RDS")
saveRDS(m.betatax100, "clean_data/biodiv_outputs/m.betatax100.RDS")
saveRDS(m.betaphy100, "clean_data/biodiv_outputs/m.betaphy100.RDS")
saveRDS(m.betafun100, "clean_data/biodiv_outputs/m.betafun100.RDS")
saveRDS(m.betapcdp100, "clean_data/biodiv_outputs/m.betapcdp100.RDS")

# or read them in...
#m.alphafun100 <- readRDS("clean_data/biodiv_outputs/m.alphafun100.RDS")
#m.alphaphy100 <- readRDS("clean_data/biodiv_outputs/m.alphaphy100.RDS")
#m.alphatax100 <- readRDS("clean_data/biodiv_outputs/m.alphatax100.RDS")
#m.betatax100 <- readRDS("clean_data/biodiv_outputs/m.betatax100.RDS")
#m.betaphy100 <- readRDS("clean_data/biodiv_outputs/m.betaphy100.RDS")
#m.betafun100 <- readRDS("clean_data/biodiv_outputs/m.betafun100.RDS")
#m.betapcdp100 <- readRDS("clean_data/biodiv_outputs/m.betapcdp100.RDS")

alpha.mam <- cbind(as.data.frame(colMeans(m.alphatax100[1:100,])), as.data.frame(colMeans(m.alphaphy100[1:100,])),
                   as.data.frame(colMeans(m.alphafun100[1:100,])))
names(alpha.mam) <- c("mammal.tax", "mammal.phy", "mammal.fun")

beta.mam <- cbind(as.data.frame(colMeans(m.betatax100[1:100,])), as.data.frame(colMeans(m.betaphy100[1:100,])),
                  as.data.frame(colMeans(m.betapcdp100[1:100,])), as.data.frame(colMeans(m.betafun100[1:100,])))
names(beta.mam) <- c("mammal.tax", "mammal.phy", "mammal.pcdp", "mammal.fun")

# calculate diversity metrics without sampling plots/dates

c(m.alphatax100.all, m.alphaphy100.all, m.alphafun100.all, m.betatax100.all, 
  m.betaphy100.all, m.betapcdp100.all, m.betafun100.all) %<-%
  .biodiv100(mammal.data, mammal.trees, m.impute, env, sample=FALSE, plots=2)

saveRDS(m.alphafun100.all, "clean_data/biodiv_outputs/m.alphafun100.all.RDS")
saveRDS(m.alphaphy100.all, "clean_data/biodiv_outputs/m.alphaphy100.all.RDS")
saveRDS(m.alphatax100.all, "clean_data/biodiv_outputs/m.alphatax100.all.RDS")
saveRDS(m.betatax100.all, "clean_data/biodiv_outputs/m.betatax100.all.RDS")
saveRDS(m.betaphy100.all, "clean_data/biodiv_outputs/m.betaphy100.all.RDS")
saveRDS(m.betafun100.all, "clean_data/biodiv_outputs/m.betafun100.all.RDS")
saveRDS(m.betapcdp100.all, "clean_data/biodiv_outputs/m.betapcdp100.all.RDS")

# or read them in...
#m.alphafun100 <- readRDS("clean_data/biodiv_outputs/m.alphafun100.all.RDS")
#m.alphaphy100 <- readRDS("clean_data/biodiv_outputs/m.alphaphy100.all.RDS")
#m.alphatax100 <- readRDS("clean_data/biodiv_outputs/m.alphatax100.all.RDS")
#m.betatax100 <- readRDS("clean_data/biodiv_outputs/m.betatax100.all.RDS")
#m.betaphy100 <- readRDS("clean_data/biodiv_outputs/m.betaphy100.all.RDS")
#m.betafun100 <- readRDS("clean_data/biodiv_outputs/m.betafun100.all.RDS")
#m.betapcdp100 <- readRDS("clean_data/biodiv_outputs/m.betapcdp100.all.RDS")

alpha.mam.all <- cbind(as.data.frame(colMeans(m.alphatax100.all[1:100,])), as.data.frame(colMeans(m.alphaphy100.all[1:100,])),
                   as.data.frame(colMeans(m.alphafun100.all[1:100,])))
names(alpha.mam.all) <- c("mammal.tax.all", "mammal.phy.all", "mammal.fun.all")

beta.mam.all <- cbind(as.data.frame(colMeans(m.betatax100.all[1:100,])), as.data.frame(colMeans(m.betaphy100.all[1:100,])),
                  as.data.frame(colMeans(m.betapcdp100.all[1:100,])), as.data.frame(colMeans(m.betafun100.all[1:100,])))
names(beta.mam.all) <- c("mammal.tax.all", "mammal.phy.all", "mammal.pcdp.all", "mammal.fun.all")

##################################
######  plant diversity    #######
##################################

plant.data <- readRDS("clean_data/PlantCommunity.RDS")
plant.trees <- readRDS("clean_data/PreppedPlantTrees.RDS")
plant.trait <- readRDS("clean_data/Clean_PlantMaxHeight.RDS")
str(env) # climate/env data from above
plant.tr <- setNames(plant.trait$height, rownames(plant.trait))
names(plant.trees) <- 1:100

# estimate phylogenetic signal and impute missing traits
c(p.phylosig100, p.impute100) %<-% .phylotrait100(p, plant.trees, plant.tr)
p.phylosig <- colMeans(p.phylosig100)
# average imputed traits
trait.mean <- as.data.frame(colMeans(p.impute100))
rownames(trait.mean) <- sort(tree.names$species)
names(trait.mean) <- "height"
saveRDS(trait.mean, "clean_data/imputedPlantTraits.rds")

# OR skip above and just read in
#plant.trait <- readRDS("clean_data/imputedPlantTraits.rds")

# 115 species from the assemblage data are missing from the tree 
# (after the congeneric.merge). These are the species:
datnames <- colnames(plant.data[4:2713])
treenames <- plant.trees[[1]]
treenames <- unique(treenames$tip.label)
setdiff(datnames, treenames)

# calculate diversity metrics with sampling plots/dates
c(p.alphatax100, p.alphaphy100, p.alphafun100, p.betatax100, 
  p.betaphy100, p.betapcdp100, p.betafun100) %<-%
  .biodiv100(plant.data, plant.trees, plant.trait, env, sample=TRUE, plot=10)

saveRDS(p.alphafun100, "clean_data/biodiv_outputs/p.alphafun100.RDS")
saveRDS(p.alphaphy100, "clean_data/biodiv_outputs/p.alphaphy100.RDS")
saveRDS(p.alphatax100, "clean_data/biodiv_outputs/p.alphatax100.RDS")
saveRDS(p.betatax100, "clean_data/biodiv_outputs/p.betatax100.RDS")
saveRDS(p.betaphy100, "clean_data/biodiv_outputs/p.betaphy100.RDS")
saveRDS(p.betafun100, "clean_data/biodiv_outputs/p.betafun100.RDS")
saveRDS(p.betapcdp100, "clean_data/biodiv_outputs/p.betapcdp100.RDS")

# or read them in...
#p.alphafun100 <- readRDS("clean_data/biodiv_outputs/p.alphafun100.RDS")
#p.alphaphy100 <- readRDS("clean_data/biodiv_outputs/p.alphaphy100.RDS")
#p.alphatax100 <- readRDS("clean_data/biodiv_outputs/p.alphatax100.RDS")
#p.betatax100 <- readRDS("clean_data/biodiv_outputs/p.betatax100.RDS")
#p.betaphy100 <- readRDS("clean_data/biodiv_outputs/p.betaphy100.RDS")
#p.betafun100 <- readRDS("clean_data/biodiv_outputs/p.betafun100.RDS")
#p.betapcdp100 <- readRDS("clean_data/biodiv_outputs/p.betapcdp100.RDS")

alpha.plant <- cbind(as.data.frame(colMeans(p.alphatax100)), as.data.frame(colMeans(p.alphaphy100)),
                     as.data.frame(colMeans(p.alphafun100)))
names(alpha.plant) <- c("plant.tax", "plant.phy", "plant.fun")

beta.plant <- cbind(as.data.frame(colMeans(p.betatax100)), as.data.frame(colMeans(p.betaphy100)),
                    as.data.frame(colMeans(p.betapcdp100)), as.data.frame(colMeans(p.betafun100)))
names(beta.plant) <- c("plant.tax", "plant.phy", "plant.pcdp", "plant.fun")

#################

# calculate diversity metrics without sampling plots/dates

c(p.alphatax100.all, p.alphaphy100.all, p.alphafun100.all, p.betatax100.all, 
  p.betaphy100.all, p.betapcdp100.all, p.betafun100.all) %<-%
  .biodiv100(plant.data, plant.trees, plant.trait, env, sample=FALSE, plot=10)

saveRDS(p.alphafun100.all, "clean_data/biodiv_outputs/p.alphafun100.all.RDS")
saveRDS(p.alphaphy100.all, "clean_data/biodiv_outputs/p.alphaphy100.all.RDS")
saveRDS(p.alphatax100.all, "clean_data/biodiv_outputs/p.alphatax100.all.RDS")
saveRDS(p.betatax100.all, "clean_data/biodiv_outputs/p.betatax100.all.RDS")
saveRDS(p.betaphy100.all, "clean_data/biodiv_outputs/p.betaphy100.all.RDS")
saveRDS(p.betafun100.all, "clean_data/biodiv_outputs/p.betafun100.all.RDS")
saveRDS(p.betapcdp100.all, "clean_data/biodiv_outputs/p.betapcdp100.all.RDS")

# or read them in...
#p.alphafun100 <- readRDS("clean_data/biodiv_outputs/p.alphafun100.all.RDS")
#p.alphaphy100 <- readRDS("clean_data/biodiv_outputs/p.alphaphy100.all.RDS")
#p.alphatax100 <- readRDS("clean_data/biodiv_outputs/p.alphatax100.all.RDS")
#p.betatax100 <- readRDS("clean_data/biodiv_outputs/p.betatax100.all.RDS")
#p.betaphy100 <- readRDS("clean_data/biodiv_outputs/p.betaphy100.all.RDS")
#p.betafun100 <- readRDS("clean_data/biodiv_outputs/p.betafun100.all.RDS")
#p.betapcdp100 <- readRDS("clean_data/biodiv_outputs/p.betapcdp100.all.RDS")

alpha.plant.all <- cbind(as.data.frame(colMeans(p.alphatax100.all)), as.data.frame(colMeans(p.alphaphy100.all)),
                     as.data.frame(colMeans(p.alphafun100.all)))
names(alpha.plant.all) <- c("plant.tax.all", "plant.phy.all", "plant.fun.all")

beta.plant.all <- cbind(as.data.frame(colMeans(p.betatax100.all)), as.data.frame(colMeans(p.betaphy100.all)),
                    as.data.frame(colMeans(p.betapcdp100.all)), as.data.frame(colMeans(p.betafun100.all)))
names(beta.plant.all) <- c("plant.tax.all", "plant.phy.all", "plant.pcdp.all", "plant.fun.all")

#####################################
###### environmental pca    #########
#####################################

### Environmental PCA
enviro <- env[,c(4:11)]
env.pca <- prcomp(enviro, scale=TRUE)
env.pc <- env.pca$x
env.pc <- env.pc[,1:2]
env.pc <- as.data.frame(env.pc)
# add elevation
env.pc$Elevation <- env$Elevation

pc.summary <- summary(env.pca)

print(xtable(pc.summary), include.rownames=TRUE, digits=c(2,2,2,2,2,2,2,2))
print(xtable(env.pca), include.rownames=TRUE)

## Environmental + Spatial Distance
pc1 <- as.numeric(dist(env.pc$PC1))
pc2 <- as.numeric(dist(env.pc$PC2))
elev <- as.numeric(dist(env$Elevation))

latlongs <- as.matrix(cbind(env$Lon,env$Lat))
distance <- distm(latlongs, fun=distHaversine)
distance[upper.tri(distance, diag=TRUE)] <- NA
dist <- as.numeric(distance)
dist <- na.omit(dist)


############################################
###### make biodiversity objects   #########
############################################

alpha.div <- cbind(alpha.plant, alpha.bird, alpha.mam, env.pc)
alpha.div$Temp <- -1*alpha.div$PC1
alpha.div <- as.data.frame(scale(alpha.div))

saveRDS(alpha.div, "clean_data/alphaDiversity.rds")

beta.mam$mammal.pcdp[beta.mam$mammal.pcdp == "NaN"] = "NA"
beta.mam$mammal.pcdp[beta.mam$mammal.pcdp == "Inf"] = "NA"
beta.mam$mammal.pcdp[beta.mam$mammal.pcdp == "-Inf"] = "NA"
beta.mam$mammal.pcdp <- as.numeric(beta.mam$mammal.pcdp)
beta.div <- cbind(beta.plant, beta.bird, beta.mam, pc1, pc2, elev, dist)
beta.div$temp <- -1*beta.div$pc1
beta.div <- as.data.frame(scale(beta.div, center=TRUE)) 

saveRDS(beta.div, "clean_data/betaDiversity.rds")

#############

# plot differences between all sites and sample averages
# ASG update these objects with the .all's
par(mfrow=c(3,3))
plant.test <- cbind(alpha.plant.all, alpha.plant)
mam.test <- cbind(alpha.mam.all, alpha.mam)
birds.test <- cbind(alpha.bird.all, alpha.bird)

# taxonomic comparison
plot(plant.test$plant.tax.all, plant.test$plant.tax, pch=20, cex=0.8, axes=FALSE, 
     xlab="all sites (tax)", ylab="BS (tax)")
axis(1)
axis(2)
abline(0,1, lwd=2, col="red")

plot(birds.test$bird.tax.all, birds.test$bird.tax, pch=20, cex=0.8, axes=FALSE, 
     xlab="all sites (tax)", ylab="BS (tax)")
axis(1)
axis(2)
abline(0,1, lwd=2, col="red")

plot(mam.test$mammal.tax.all, mam.test$mammal.tax, pch=20, cex=0.8, axes=FALSE, 
     xlab="all sites (tax)", ylab="BS (tax)")
axis(1)
axis(2)
abline(0,1, lwd=2, col="red")


#### phylogenetic comparison
plot(plant.test$plant.phy.all, plant.test$plant.phy, pch=20, cex=0.8, axes=FALSE, 
     xlab="all sites (phy)", ylab="BS (phy)")
axis(1)
axis(2)
abline(0,1, lwd=2, col="red")

plot(birds.test$bird.phy.all, birds.test$bird.phy, pch=20, cex=0.8, axes=FALSE, 
     xlab="all sites (phy)", ylab="BS (phy)")
axis(1)
axis(2)
abline(0,1, lwd=2, col="red")

plot(mam.test$mammal.phy.all, mam.test$mammal.phy, pch=20, cex=0.8, axes=FALSE, 
     xlab="all sites (phy)", ylab="BS (phy)")
axis(1)
axis(2)
abline(0,1, lwd=2, col="red")

#### functional comparison
plot(plant.test$plant.fun.all, plant.test$plant.fun, pch=20, cex=0.8, axes=FALSE, 
     xlab="all sites (fun)", ylab="BS (fun)")
axis(1)
axis(2)
abline(0,1, lwd=2, col="red")

plot(birds.test$bird.fun.all, birds.test$bird.fun, pch=20, cex=0.8, axes=FALSE, 
     xlab="all sites (fun)", ylab="BS (fun)")
axis(1)
axis(2)
abline(0,1, lwd=2, col="red")

plot(mam.test$mammal.fun.all, mam.test$mammal.fun, pch=20, cex=0.8, axes=FALSE, 
     xlab="all sites (fun)", ylab="BS (fun)")
axis(1)
axis(2)
abline(0,1, lwd=2, col="red")
