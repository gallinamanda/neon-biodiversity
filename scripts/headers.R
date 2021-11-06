change.rownames <- function(data, column){
  rownames(data) <- data[,column]
  data <- data[,!names(data) %in% column,drop=FALSE]
  return(data)
}

env.average <- function(data, coords){
  val <- raster:::extract(x=data, y=coords, df=TRUE)
  recent <- val[,c(1, 90:119)]
  average <- as.data.frame(rowMeans(recent))
  return(average)
}

# make a distance matrix
convert.dist <- function(data){
  mat <- diag(38)
  mat[upper.tri(mat, diag=FALSE)] <- data
  #mat[lower.tri(mat)] <- NA
  mat <- t(mat)
  mat <- as.dist(mat, diag = FALSE, upper = FALSE)
  return(mat)
}

#generate comparative community objects
.format.c.data <- function(data, tree, traits, env, p.a=TRUE){
  comm <- picante::sample2matrix(data[,c("site","value","species")])
  comm <- data.matrix(comm, rownames.force=NA)
  
  # Convert to presence/absence
  if(p.a){
    for(i in seq_len(ncol(comm)))
      comm[comm[,i] != 0,i] <- 1
  }
  
  # Make longform data from comparative.comm
  c.data <- comparative.comm(tree, comm)
  l.data <- as.data.frame(c.data)
  
  # Add traits and environmental data
  c.data.tr <- comparative.comm(tree, comm, traits)
  l.data.tr <- as.data.frame(c.data.tr)
  c.data.full <- comparative.comm(tree, comm, traits, env)
  l.data.full <- as.data.frame(c.data.full)
  
  return(list(c.comm=c.data, l.data=l.data, l.data.full=l.data.full, c.data.full=c.data.full))
}

# NACDB functions (and utility) for downloading NEON data

#utility:

.df.melt <- function(species, site.id, value,
                     study.metadata=data.frame(units=NA, other=NA),
                     site.metadata=data.frame(id=NA,year=NA,name=NA,lat=NA,long=NA,address=NA,area=NA,other=NA),
                     species.metadata=data.frame(species=NA, taxonomy=NA, other=NA)){
  #######################
  # Argument handling ###
  #######################
  if(!is.numeric(value))
    stop("'value' is not numeric")
  if(any(is.na(value)))
    stop("No NAs in 'value'")
  if(any(is.na(species)))
    stop("No NAs in 'species'")
  if(any(is.na(site.id)))
    stop("No NAs in 'site.id'")
  species <- as.character(species)
  site.id <- as.character(site.id)
  
  ######################
  # Meta-data ##########
  ######################
  .create.other <- function(metadata, columns){
    if(!all(names(metadata) %in% columns)){
      other <- metadata[,!names(metadata) %in% columns, drop=FALSE]
      metadata <- metadata[,names(metadata) %in% columns, drop=FALSE]        
      other <- sapply(seq_along(names(other)), function(y) paste(names(other)[y],other[,y],sep=":"))
      if(ncol(metadata) == 1)
        other <- paste(other, collapse=";") else other <- apply(other, 1, paste, collapse=";")
      metadata$other <- other
    } else {
      metadata$other <- NA
    }
    return(metadata)
  }
  # Study
  if(nrow(study.metadata) > 1)
    stop("Only one row of meta-data per study")
  if(!all("units" %in% names(study.metadata)))
    stop("Incorrectly formatted study meta-data")
  if(is.na(study.metadata$units))
    stop("Study must have units of measurement")
  study.metadata <- .create.other(study.metadata, "units")
  # Site
  if(!all(c("id","year","name","lat","long","address","area") %in% names(site.metadata)))
    stop("Incorrectly formatted site meta-data")
  if(length(intersect(unique(site.id), site.metadata$id)) != nrow(site.metadata))
    stop("Site meta-data must contain information about all sites")
  if(length(intersect(site.metadata$id,unique(site.id))) != nrow(site.metadata))
    stop("Site meta-data must only contain information about present sites")
  site.metadata <- .create.other(site.metadata, c("id","year","name","lat","long","address","area"))
  # Species
  if(!all(c("species","taxonomy") %in% names(species.metadata)))
    stop("Incorrectly formatted species meta-data")
  if(length(intersect(unique(species), species.metadata$species)) != nrow(species.metadata))
    stop("Species meta-data must contain information about all species")
  if(length(intersect(species.metadata$species,unique(species))) != nrow(species.metadata))
    stop("Species meta-data must only contain information about present species")
  species.metadata <- .create.other(species.metadata, c("species","taxonomy"))
  
  ######################
  # Format and return ##
  ######################
  # Reformat data
  output <- list(
    data=data.frame(site.id, species, value),
    spp.metadata=species.metadata,
    site.metadata=site.metadata,
    study.metadata=study.metadata
  )
  for(i in seq_along(output))
    for(j in seq_len(ncol(output[[i]])))
      if(is.factor(output[[i]][,j]))
        output[[i]][,j] <- as.character(output[[i]][,j])
  class(output) <- "nacdb"
  return(output)
}


# Mammals
#' @importFrom neonUtilities loadByProduct
#' @export
.neon.2018a <- function(...){
  dat <- loadByProduct("DP1.10072.001", site ="all", startdate = "2017-01", enddate = "2017-12")
  trapnight <- dat$mam_pertrapnight
  trapnight$id <- paste0(trapnight$plotID, "_", trapnight$collectDate)
  data <- trapnight[,c("id", "plotID", "collectDate", "decimalLatitude", "decimalLongitude", "scientificName"),]
  data <- data[which(data$scientificName!=""),]
  names(data) <- c("id", "name", "year", "lat", "long", "species")
  temp <- with(data, paste(id, species))
  temp.counts <- table(temp)
  data$value <- temp.counts[temp]
  data$value <- as.numeric(data$value)
  data <- unique(data)
  data <- na.omit(data)
  site.df <- data[, c("id", "name", "year", "lat", "long")]
  site.df <- unique(site.df)
  site.df$address <- NA; site.df$area <- "10 sherman traps"
  return(.df.melt(data$species, data$id, data$value, data.frame(units="g"), site.df, data.frame(species=unique(data$species),taxonomy=NA)))
}

# Plants
#' @export
.neon.2018c <- function(...){
  dat <- loadByProduct("DP1.10058.001", site ="all", startdate = "2017-01", enddate = "2017-12")
  vst <- dat$div_10m2Data100m2Data
  vst$date <- vst$endDate
  vst$id <- paste0(vst$plotID, "_", vst$date)
  data <- vst[,c("id", "plotID", "date", "scientificName"),]
  names(data) <- c("id", "name", "year", "species")
  data$species <- str_extract(data$species, "[^ ]+ [^ ]+")
  temp <- with(data, paste(id, species))
  temp.counts <- table(temp)
  data$value <- temp.counts[temp]
  data$value <- as.numeric(data$value)
  data <- unique(data)
  data <- na.omit(data)
  site.df <- data[, c("id", "name", "year")]
  site.df <- unique(site.df)
  site.df$lat <- NA; site.df$long <- NA
  site.df$address <- NA; site.df$area <- "1m2"
  return(.df.melt(data$species, data$id, data$value, data.frame(units="#"), site.df, data.frame(species=unique(data$species),taxonomy=NA)))
}

# Birds
#' @export
.neon.2018d <- function(...){
  dat <- loadByProduct("DP1.10003.001", site ="all", startdate = "2017-01", enddate = "2017-12")
  vst <- dat$brd_countdata
  vst$date <- substr(vst$startDate,0,10)
  vst$id <- paste0(vst$plotID, "_", vst$date)
  data <- vst[,c("id", "plotID", "date", "scientificName", "clusterSize"),]
  names(data) <- c("id", "name", "year", "species", "value")
  data <- aggregate(data$value, by=list(data$id, data$name, data$year, data$species), FUN=sum, na.rm=TRUE)
  names(data) <- c("id", "name", "year", "species", "value")
  site.df <- data[, c("id", "name", "year")]
  site.df <- unique(site.df)
  site.df$lat <- NA; site.df$long <- NA
  site.df$address <- NA; site.df$area <- "?"
  return(.df.melt(data$species, data$id, data$value, data.frame(units="#"), site.df, data.frame(species=unique(data$species),taxonomy=NA)))
}


# check phylogenetic signal, then impute across trees and average

.phylotrait100 <- function(data, trees, traits){
  tree.list <- unique(names(trees))
  phylosig100 <- matrix(NA, nrow=length(trees), ncol=2)
  impute100 <- list()
  for(k in 1:length(tree.list)){
  tree <- trees[[k]]
  c(c.dat, l.dat, l.dat.env, c.dat.env) %<-% 
    .format.c.data(data, tree, trait=NULL, env=NULL) 
  post.tree <- c.dat$phy
  sig <- phylosig(post.tree, traits, method="lambda", test=TRUE)
  phylosig100[k,] <- cbind(sig$lambda, sig$P)
  
  tree.names <- as.data.frame(post.tree$tip.label); names(tree.names) <- "species"
  trait <- as.data.frame(traits); trait$species <- rownames(trait)
  new.trait <- merge(trait, tree.names, by="species", all.y=TRUE)
  names(new.trait) <- c("species", "value")
  phy <- phylopars(new.trait, post.tree, model="BM")
  recon <- as.data.frame(phy$anc_recon)
  recon$species <- rownames(recon)
  rownames(recon) <- NULL
  recon <- merge(recon, tree.names, by="species")
  imputed <- as.data.frame(t(recon$value))
  impute.names <- recon$species
  impute100 <- rbind(impute100, imputed)
  }
  names(impute100) <- impute.names
  return(list(phylosig100=phylosig100, impute100=impute100))
}


# diversity calculations-- sampling from plots/dates is optional
.biodiv100 <- function(data, trees, traits, env, sample=TRUE, plots){
  tree.list <- unique(names(trees))
  alphatax100 <- matrix(NA, nrow=100, ncol=38)
  alphaphy100 <- matrix(NA, nrow=100, ncol=38)
  alphafun100 <- matrix(NA, nrow=100, ncol=38)
  betatax100 <- matrix(NA, nrow=100, ncol=703)
  betaphy100 <- matrix(NA, nrow=100, ncol=703)
  betapcdp100 <- matrix(NA, nrow=100, ncol=703)
  betafun100 <- matrix(NA, nrow=100, ncol=703)
  
  for(k in 1:length(tree.list)){
    tree <- trees[tree.list[[k]]]
    tree <- tree[[1]]
    #for loop by site
    sitelist<-unique(as.character(data$site))
    all.sites <- NULL
    all.sites <- list()
    for(j in 1:length(sitelist)){
      
      #subset to one site
      single.site <- subset(data, site==sitelist[j])
      
      if(sample){
      #randomly select one unique row for each plot, so that only one row for each plot remains
      single.site <- stratified(single.site, group= "plot", size=1)
      # randomly select minimum number of plots sampled
      single.site <- stratified(single.site, group= "site", size=plots)
      }
      
      single.site <- single.site[,4:ncol(single.site)]
      one.site <- as.data.frame(t(colSums(single.site)))
      site <- sitelist[j]
      
      # add the name of the site as a column
      one.site <- as.data.frame(cbind(as.data.frame(site), one.site))
      
      # do this for all sites and make that into a data frame
      all.sites<-rbind(all.sites, one.site)
    }
    
    # first make all.sites into a DF of site, species, value
    all.sites <- all.sites %>% mutate_if(is.numeric, ~1 * (. > 0)) # make into occurrence by site
    sites_long <- gather(all.sites, species, value, 
                         colnames(all.sites[2]):colnames(all.sites[ncol(all.sites)]), factor_key=TRUE)
    
    # first, have to limit traits to species present at *some* site, or else an error:
    # "At least one species does not occur in any community (zero total abundance 
    # across all communities)"
    sp.dat <- all.sites[,c(2:ncol(all.sites))]
    z <- (colSums(sp.dat, na.rm=T) != 0)
    nonzerospp <- sp.dat[,z]
    nonzerospp <- cbind(all.sites$site, nonzerospp)
    nonzero_long <- gather(nonzerospp, species, value, 
                           colnames(nonzerospp[2]):colnames(nonzerospp[ncol(nonzerospp)]), factor_key=TRUE)
    names(nonzero_long)[1] <- "site"
    
    # make a comparative comm
    c(c.dat, l.dat, l.dat.env, c.dat.env) %<-% 
      .format.c.data(nonzero_long, tree, traits, env)
    
    # use the comparative comm to estimate the following diversity metrics--
    comm <- c.dat$comm
    phy <- c.dat$phy
    
    # alpha diversity: taxonomic (species richness)
    alphatax100[k,] <- as.numeric(rowSums(comm))
    
    # alpha diversity: phylogenetic (ses.mntd)
    phylo <- cophenetic(phy)
    alphaphy <- ses.mntd(comm, phylo)
    alphaphy100[k,] <- as.numeric(alphaphy$mntd.obs.z)
    
    # beta diversity: taxonomic (sorensen index)
    betatax100[k,] <- as.numeric(vegdist(comm, method="bray", binary=TRUE))
    
    #pez dissimilarity:
    pez <- pez.dissimilarity(c.dat)
    # beta diversity: phylogenetic (phylosorensen index)
    betaphy100[k,] <- as.numeric(pez$phylosor)
    
    # beta diversity: phylogenetic (PCDp)
    betapcdp100[k,] <- as.numeric(pez$pcd$PCDp)
    
    # alpha diversity: functional (dispersion)
    alphafun <- dbFD(c.dat.env$data, c.dat.env$comm, calc.CWM=TRUE)
    alphafun100[k,] <- as.numeric(alphafun$FDis)
    
    # beta diversity: functional (dissimilarity)
    sptrait <- traits.dist(c.dat.env)
    betafun <- pez.dissimilarity(c.dat.env, ext.dist = sptrait)
    betafun100[k,] <- as.numeric(betafun$comdist)
  }
  
  return(list(alphatax100=alphatax100, alphaphy100=alphaphy100, 
              alphafun100=alphafun100, betatax100=betatax100,
              betaphy100=betaphy100, betapcdp100=betapcdp100, 
              betafun100=betafun100))
}

coefs <- function(model){
  coefs <- as.data.frame(model$coefficients)
  coefs <- coefs[c(2:6),c(1,4)]
  names(coefs) <- c("est", "p")
  est <- coefs$est
  p <- coefs$p
  #ns <- within(coefs, est[p <=0.05] <- 0)
  #ns <- within(ns, est[p >=0.05] <- 0.2)
  #ns <- ns$est
  return(list(est=est, p=p))
}

b.coefs <- function(model){
  coefs <- as.data.frame(model$coefficients)
  coefs <- coefs[c(2:7),c(1,4)]
  names(coefs) <- c("est", "p")
  est <- coefs$est
  p <- coefs$p
  return(list(est=est, p=p))
}
