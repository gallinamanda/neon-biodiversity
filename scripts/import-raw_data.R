# Load (and install, if needed) packages and convenience functions
source("scripts/setup.R")
source("scripts/headers.R")

###############################
########## Plant Data #########
###############################

###### Download and save plant inventories from NEON
plant <- .neon.2018c()
saveRDS(plant, "raw_data/NEONplant_raw_data.rds")

###### Plant Traits from BIEN
# prep the species list
plant.names <- as.data.frame(unique(plant$data$species))
names(plant.names) <- "species"
# remove all names that are genus only (sp. and spp. removed)
keep <- as.data.frame(plant.names[!grepl("\\.", plant.names$species),])
names(keep) <- "species"
# remove unknowns
keep <- keep[(!(keep$species=="Unknown plant") & !(keep$species=="Unknown softwood") 
                & !(keep$species=="Unknown hardwood")),]
keep <- as.character(keep)
saveRDS(keep, "clean_data/CleanPlantNames.RDS")

# download plant height data for species list
plant.height <- c("whole plant height")
plant.height <- BIEN_trait_traitbyspecies(species = keep, trait = plant.height, all.taxonomy=FALSE)
plant.maxheight <- data.frame(plant.height$scrubbed_species_binomial, as.numeric(plant.height$trait_value))
names(plant.maxheight) <- c("species", "plant_height_m")
# calculate maximum plant height
plant.maxheight <- aggregate(. ~ species, data=plant.maxheight, max)
plant.maxheight$species <- gsub(" ", "_", plant.maxheight$species)
# remove plant height = 0 from data set
plant.maxheight <- plant.maxheight[plant.maxheight$plant_height_m != 0,]
# calculate log(max plant height)
plant.maxheight$plant_height_m <- as.numeric(log(plant.maxheight$plant_height_m))
names(plant.maxheight) <- c("species", "log.maxheight")
plant.traitspp <- as.data.frame(unique(plant.maxheight$species)) #977 species represented

# save out to raw_data folder
saveRDS(plant.maxheight, "raw_data/PlantTraits.rds")

###############################
########## Mammal Data ########
###############################

###### Download and save mammal inventories from NEON
mammal <- .neon.2018a()
saveRDS(mammal, "raw_data/NEONMammal_raw_data.rds")

###### Mammal Traits from EltonTraits
mam.trait <- fread('http://www.esapubs.org/archive/ecol/E095/178/MamFuncDat.txt')
mam.trait <- data.frame(mam.trait$Scientific, mam.trait$`BodyMass-Value`)
names(mam.trait) <- c("species", "log.bodymass")
mam.trait <- na.omit(mam.trait)
# calculate log(mammal body mass)
mam.trait$log.bodymass <- as.numeric(log(mam.trait$log.bodymass))
mam.trait$species <- gsub(" ", "_", mam.trait$species)

# save out to raw_data folder
saveRDS(mam.trait, "raw_data/MammalTraits.rds")

###############################
########## Bird Data #########
###############################

###### Download and save bird inventories from NEON
bird <- .neon.2018d()
saveRDS(bird, "raw_data/NEONbird_raw_data.rds")

###### Bird Traits from EltonTraits
bird.trait <- fread('http://www.esapubs.org/archive/ecol/E095/178/BirdFuncDat.txt')
bird.trait <- data.frame(bird.trait$Scientific, bird.trait$`BodyMass-Value`)
names(bird.trait) <- c("species", "log.bodymass")
bird.trait <- na.omit(bird.trait)
# calculate log(bird body mass)
bird.trait$log.bodymass <- as.numeric(log(bird.trait$log.bodymass))
bird.trait$species <- gsub(" ", "_", bird.trait$species)

# save out to raw_data folder
saveRDS(bird.trait, "raw_data/BirdTraits.rds")

###############################
## Environment/Site Data ######
###############################

#### Read in environmental data
#CSV of basic site info from NEON website
neon.field <- read.csv("raw_data/field-sites.csv")
field.sites <- neon.field[,c("Site.ID", "Domain.Number", "State", "Latitude", "Longitude", 
                             "Elevation")]
field.sites$Elevation <- as.numeric(gsub("[^0-9]+", "",  field.sites$Elevation))
names(field.sites) <- c("site", "Domain", "State", "Lat", "Lon", "Elevation")

# CRU Climate Data
meantemp <- readRDS("raw_data/cru-tmp.rds")
maxtemp <- readRDS("raw_data/cru-tmx.rds")
mintemp <- readRDS("raw_data/cru-tmn.rds")
precip <- readRDS("raw_data/cru-pre.rds")
vapor <- readRDS("raw_data/cru-vap.rds")
potential <- readRDS("raw_data/cru-pet.rds")
cloud <- readRDS("raw_data/cru-cld.rds")
frost <- readRDS("raw_data/cru-frs.rds")
  
# Process and collate environmental data
sitedata <- field.sites
sitedata <- sitedata[,c("site", "Lon", "Lat")]
names(sitedata) <- c("site", "lon", "lat")
neon.coords <- data.frame(lon=sitedata[,2], lat=sitedata[,3])
coordinates(neon.coords) <- c("lon", "lat")

# Add meteorological variables
env <- cbind(sitedata$site, env.average(meantemp, neon.coords), env.average(mintemp, neon.coords), 
             env.average(maxtemp, neon.coords), env.average(precip, neon.coords), env.average(vapor, neon.coords), 
             env.average(potential, neon.coords), env.average(cloud, neon.coords), env.average(frost, neon.coords)) 

names(env) <- c("site", "tmp", "tmn", "tmx", "pre", "vap", "pet", "cld", "frs")
climate.env <- merge(field.sites, env, by="site")

# save out to raw_data folder
saveRDS(climate.env, "raw_data/climate_data.RDS")
