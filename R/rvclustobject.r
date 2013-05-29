#' Constructor for rvclustobject
#'
#' ---------------------------------------------------------------------- '#\cr
#' RVCLUST functionality is modularized for easy extension development. \cr
#' Each module is standardized to accept an rvclustobject, which is \cr
#' initialized with variant, observation, covariate, and outcome data \cr
#' files. \cr
#' \cr
#' Variant annotations may be applied by \code{\link{annotate}}. \cr
#' \cr
#' RVCLUST was originally designed for the analysis of rare genetic \cr
#' variants. In line with its original intent, PEDMAP files may be \cr
#' used as input as an alternative to the standardized input format. \cr
#' This functionality is dependent on PLINK. \cr
#' ---------------------------------------------------------------------- '#\cr
#' 
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param observation.file observation data (or PED)
#' @param variant.file variant data (or MAP)
#' @param cov.file covariate data
#' @param outcome.file outcome data
#' @return rvclustobject with fields\cr
#'  \tabular{ll}{
#'  variants\cr
#'  observations\cr
#'  covariates\cr
#'  outcomes\cr
#'  clusters\cr
#'  clusterinfo\cr
#'  collapsed\cr
#'  }
rvclustobject <- function(observation.file=NA,variant.file=NA,cov.file=NA,outcome.file=NA,max.freq=100.0) {
  if (is.na(observation.file)) {
    rv <- load_sample()}

  else if (tolower(unlist(strsplit(basename(observation.file),split="\\."))[2]) == "ped") {
    rv <- load_pedmap(observation.file,variant.file,cov.file,outcome.file,max.freq)}
  
  else {
    rv <- load_default(observation.file,variant.file,cov.file,outcome.file,max.freq)}
  return(rv)
}

show <- function(rv) {UseMethod("show",rv)}
show.rvclustobject <- function(rv) {print(c("data","variants","clusters","clusterinfo"))}

print <- function(rv) {UseMethod("print",rv)}
print.rvclustobject <- function(rv) {print(c("data","variants","clusters","clusterinfo"))}

summary <- function(rv) {UseMethod("summary",rv)}
summary.rvclustobject <- function(rv) {
  print("rvclust object:")
  print("Data:")
  print(summary(rv$observations))
  print(summary(rv$variants))
  print(summary(rv$covariates))
  print(summary(rv$outcomes))
  print("Variants:")
  print(summary(rv$variants))
  print("Clusters:")
  print(summary(rv$clusters))
  print("Cluster Info:")
  print(summary(rv$clusterinfo))
}

load_sample <- function() {
  data(map.dat)
  names(map.dat) <- c('CHR','SNP','DIST','POS')
  data(raw.dat)
  names(raw.dat) <- sapply(names(raw.dat),function(x){x <- strsplit(x,'_')[[1]][1]; sub("[[punct:]]",".",x)})
  rv <- list(
            "observations" =raw.dat,
            "covariates"   =NA,
            "outcomes"     =NA,
            "variants"     =rv.dat,
            "clusters"     =NA,
            "clusterinfo"  =NA,
            "collapsed"=NA)
  class(rv) <- "rvclustobject"
  return(rv)
}

load_default <- function(observation.file=NA,variant.file=NA,cov.file=NA,outcome.file=NA,max.freq=100.0) {
  observations <- read.table(observation.file,sep='',header=TRUE)
  variants     <- read.table(variant.file,sep='',header=TRUE)
  if (!is.na(cov.file)) {
    covariates <- read.table(cov.file,sep='',header=TRUE)}
  else{covariates=NA}
  if (!is.na(outcome.file)) {
    outcomes   <- read.table(outcome.file,sep='',header=TRUE)}
  else{outcomes=NA}
  
  if (max.freq < 100.00) {
    writeLines("RVCLUST does not yet support frequency calculations for non-genetic data.")}

  rv <- list(
            "observations" =observations,
            "variants"     =variants,
            "covariates"   =covariates,
            "outcomes"     =outcomes,
            "clusters"     =NA,
            "clusterinfo"  =NA,
            "collapsed.dat"=NA)
  class(rv) <- "rvclustobject"
  return(rv)
}

load_pedmap <- function(ped.file=NA,map.file=NA,cov.file=NA,phen.file=NA,max.freq=100.0) {
  pedmap.path  <- dirname(ped.file)
  pedmap.fname <- unlist(strsplit(basename(ped.file),split="\\."))[1]
  raw.file  <- paste(pedmap.path,"/",pedmap.fname,".raw",sep='')

  if (!file.exists(raw.file)) {
    system(paste("plink --noweb --file",paste(pedmap.path,"/",pedmap.fname,sep=''),"--recodeA --out",paste(pedmap.path,"/",pedmap.fname,sep='')))}
    
  map.dat  <- read.table(map.file,sep='',header=FALSE,col.names=c('CHR','SNP','DIST','POS'))
  map.dat$SNP <- sapply(map.dat$SNP,function(x){sub("[[:punct:]]",".",x)})
  raw.dat  <- read.table(raw.file,sep='',header=TRUE)
  names(raw.dat) <- sapply(names(raw.dat),function(x){x <- strsplit(x,'_')[[1]][1]; sub("[[punct:]]",".",x)})
  cov.dat  <- NA
  if (!is.na(cov.file)) {
    cov.dat <- read.table(cov.file,sep='',header=TRUE)}
  phen.dat <- NA
  if (!is.na(phen.file)) {
    phen.dat <- read.table(phen.file,sep='',header=TRUE)}

  rv.dat <- rare.vars(pedmap.path,pedmap.fname,map.dat,max.freq)

  rv <- list(
            "observations" =raw.dat,
            "covariates"   =cov.dat,
            "outcomes"     =phen.dat,
            "variants"     =rv.dat,
            "clusters"     =NA,
            "clusterinfo"  =NA,
            "collapsed.dat"=NA)
  class(rv) <- "rvclustobject"
  return(rv)
}

# --------------------------------------------------------------------------- #
# Deprecated methods below are required for PEDMAP processing.
# --------------------------------------------------------------------------- #

rare.vars <- function(pedmap.path,pedmap.fname,map.dat,max.freq) {
  # Identify the major/minor allele in all SNPs

  # If running an example
  if (is.na(pedmap.path) & is.na(pedmap.fname)) {
    data(rv.dat)
  }
  
  else {
    snp.dat <- allele.frequency(pedmap.path,pedmap.fname)
    snp.dat$SNP <- sapply(snp.dat$SNP,function(x){sub("[[:punct:]]",".",x)})

    # Separate the common and rare variants
    rv.dat <- rare.snps(snp.dat,max.freq)
    
    # Use the map file to add position information
    rv.dat <- add.pos(rv.dat,map.dat)
  }

  return(rv.dat)
}

allele.frequency <- function(pedmap.path,pedmap.fname) {
  # Given a PED/MAP path and basename, use PLINK to
  # return the rsIDs with frequency > 0 and < 0.05
  
  # Build the full file
  fin  <- paste(pedmap.path,'/',pedmap.fname,sep='')
  
  # Call PLINK to generate frequency file
  system(paste("plink --noweb --file",fin,"--freq --allow-no-sex",sep=' '))
  
  # Give PLINK time to write the file
  Sys.sleep(2)
  
  # Read the frequency file
  freq.dat <- read.table("plink.frq",sep='',header=TRUE,col.names=c("CHR","SNP","MA","MajA","MAF","NCHROBS"))
    
  # Subset the data to include only SNP, Minor Allele, and Frequency info
  freq.dat <- subset(freq.dat,select=c("SNP","MA","MAF"))
}

rare.snps <- function(freq.dat,max.freq) {
  # Removes all common variants and returns the rare subset
  rare.snp.dat <- subset(freq.dat,MAF>0 & MAF<max.freq)
}

add.pos <- function(rv.dat,map.dat) {
  # Uses the MAP file to add position info to the rare variant data frame
  
  col.names <- append(names(rv.dat),c("CHR","POS"))
  rv.dat <- merge(rv.dat,map.dat,by="SNP")
  rv.dat <- subset(rv.dat,select=col.names)
}