###################################
# Installing, and then loading, ###
# all required packages         ###
###################################
# Note: there is an intention of doing these installs in ~vaguely an
# order of their own dependencies (e.g., install ape first because
# everything else depends on it)

# Phylogenetics
if(!require(ape)){
    install.packages("ape")
    library(ape)
}
if(!require(phytools)){
    install.packages("phytools")
    library(phytools)
}
if(!require(picante)){
    install.packages("picante")
    library(picante)
}
if(!require(pez)){
    install.packages("pez")
    library(pez)
}
if(!require(Rphylip)){
    install.packages("Rphylip")
    library(Rphylip)
}

if(!require(vegan)){
  install.packages("vegan")
  library(vegan)
}

# Data manipulation
if(!require(dplyr)){
    install.packages("dplyr")
    library(dplyr)
}
if(!require(stringr)){
    install.packages("stringr")
    library(stringr)
}
if(!require(doBy)){
    install.packages("doBy")
    library(doBy)
}
if(!require(data.table)){
    install.packages("data.table")
    library(data.table)
}
if(!require(tidyr)){
    install.packages("tidyr")
    library(tidyr)
}
if(!require(zeallot)){
    install.packages("zeallot")
    library(zeallot)
}

# Misc
if(!require(bit64)){
    install.packages("bit64")
    library(bit64)
}
if(!require(raster)){
    install.packages("raster")
    library(raster)
}
if(!require(BIEN)){
    install.packages("BIEN")
    library(BIEN)
}
if(!require(lme4)){
    install.packages("lme4")
    library(lme4)
}
if(!require(xtable)){
    install.packages("xtable")
    library(xtable)
}
if(!require(sjstats)){
    install.packages("sjstats")
    library(sjstats)
}
if(!require(openintro)){
    install.packages("openintro")
    library(openintro)
}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(gridExtra)){
  install.packages("gridExtra")
  library(gridExtra)
}
if(!require(lattice)){
  install.packages("lattice")
  library(lattice)
}
if(!require(corrplot)){
  install.packages("corrplot")
  library(corrplot)
}
if(!require(maptools)){
  install.packages("maptools")
  library(maptools)
}
if(!require(rgdal)){
  install.packages("rgdal")
  library(rgdal)
}
if(!require(raster)){
  install.packages("raster")
  library(raster)
}
if(!require(rgeos)){
  install.packages("rgeos")
  library(rgeos)
}
if(!require(splitstackshape)){
  install.packages("splitstackshape")
  library(splitstackshape)
}



# GitHub (and devtools dependency)
if(!require(devtools)){
    install.packages("devtools")
    library(devtools)
}
if(!require(suppdata)){
    install_github("ropensci/suppdata")
    library(suppdata)
}
if(!require(neonUtilities)){
    install.packages("neonUtilities")
    library(neonUtilities)
}
if(!require(FD)){
  install.packages("FD")
  library(FD)
}
if(!require(geosphere)){
  install.packages("geosphere")
  library(geosphere)
}

if(!require(quantreg)){
  install.packages("quantreg")
  library(quantreg)
}

if(!require(DataCombine)){
  install.packages("DataCombine")
  library(DataCombine)
}

