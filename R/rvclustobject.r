#' Constructor for rvclustobject
#'
#' Constructor for rvclustobject. Loads ped/map data from the\cr
#' specified path and filename. If a covariate path and\cr
#' filename are specified, the covariates will be loaded as well.\cr
#' Burden testing can be enabled by setting the burden flag to true.\cr
#' A threshold can also be set such that only clusters with fitness\cr
#' exceeding the threshold are tested for association.\cr
#' \cr
#' The rvclust package contains wrappers for various clustering\cr
#' methods and statistical tests. See package documentation for more\cr
#' information on these wrappers. Simply pass the rvclustobject to any\cr
#' method with an rvclust wrapper and an updated rvclustobject will be returned.
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param pedmap.path path to pedmap files
#' @param pedmap.fname basename for pedmap files (no extension)
#' @param cov.path path to covariates file [optional]
#' @param cov.fname filename for covariates file (with extension) [optional]
#' @param annotations list of desired annotations [default ALL]
#' @param burden boolean indicating whether to burden test [default FALSE]
#' @param min.fit minimum cluster fitness to test (0..1) [default 0.0]
#' @return rvclustobject containing\cr
#'  \tabular{ll}{
#'  data: \tab ped/map/covariate data\cr
#'  variants: \tab variant data.frame\cr
#'  clusters: \tab empty list of clusters\cr
#'  clusterinfo: \tab empty list of cluster info\cr
#'  }
#' @note Data must be in PEDMAP format
#' @note Covariates - One column for subject; One column per covariate
#' @examples
#' rvclustobject(NA,NA,annotations=c("CHROMATIN"))
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

# Deprecated. Required for native PEDMAP load.
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

# Deprecated. Required for native PEDMAP load.
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

# Deprecated. Required for native PEDMAP load.
add.pos <- function(rv.dat,map.dat) {
  # Uses the MAP file to add position info to the rare variant data frame
  
  col.names <- append(names(rv.dat),c("CHR","POS"))
  rv.dat <- merge(rv.dat,map.dat,by="SNP")
  rv.dat <- subset(rv.dat,select=col.names)
}

# To be removed
load.data <- function(pedmap.path,pedmap.fname,cov.path,cov.fname,phen.path,phen.fname) {

  # Sample run
  if (is.na(pedmap.path) & is.na(pedmap.fname)) {
    data(map.dat)
    names(map.dat) <- c('CHR','SNP','DIST','POS')
    data(raw.dat)
    names(raw.dat) <- sapply(names(raw.dat),function(x){x <- strsplit(x,'_')[[1]][1]; sub("[[punct:]]",".",x)})
    cov.dat <- NA
    phen.dat <- NA
  }  
  else {
    # Initialize all file names
    ped.file  <- paste(pedmap.path,"/",pedmap.fname,".ped",sep='')
    map.file  <- paste(pedmap.path,"/",pedmap.fname,".map",sep='')
    raw.file  <- paste(pedmap.path,"/",pedmap.fname,".raw",sep='')
    cov.file  <- NA
    phen.file <- NA
    if (!is.na(cov.path) & !is.na(cov.fname)) {
      cov.file <- paste(cov.path,"/",cov.fname,sep='')
    }
    if (!is.na(phen.path) & !is.na(phen.fname)) {
      phen.file <- paste(phen.path,"/",phen.fname,sep='')
    }
    
    # If the RAW file is missing, use PLINK to create it
    if (!file.exists(raw.file)) {
      system(paste("plink --noweb --file",paste(pedmap.path,"/",pedmap.fname,sep=''),"--recodeA --out",paste(pedmap.path,"/",pedmap.fname,sep='')))
    }
    
    # Load the data frames
    map.dat  <- read.table(map.file,sep='',header=FALSE,col.names=c('CHR','SNP','DIST','POS'))
    map.dat$SNP <- sapply(map.dat$SNP,function(x){sub("[[:punct:]]",".",x)})
    raw.dat  <- read.table(raw.file,sep='',header=TRUE)
    names(raw.dat) <- sapply(names(raw.dat),function(x){x <- strsplit(x,'_')[[1]][1]; sub("[[punct:]]",".",x)})
    cov.dat  <- NA
    if (!is.na(cov.file)) {
      cov.dat <- read.table(cov.file,sep='',header=TRUE)
    }
    phen.dat <- NA
    if (!is.na(phen.file)) {
      phen.dat <- read.table(phen.file,sep='',header=TRUE)
    }
  }
  return(list(map.dat,raw.dat,cov.dat,phen.dat)) 
}
