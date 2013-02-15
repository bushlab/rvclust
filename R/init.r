#' Initializes data and covariates for rare variant analysis
#'
#' Utility to load data, phenotype, and an optional covariates\cr
#' file into data frames formatted for rare variant clustering\cr
#' and analysis with rvclust.
#'
#' @export
#' @param data.path data directory
#' @param data.fname data file name
#' @param cov.path covariate file directory
#' @param cov.fname covariate file name [optional]
#' @return rvclust object containing\cr
#'  \tabular{ll}{
#'  rv.dat: \tab rare variant data frame\cr
#'  raw.dat: \tab genotype/phenotype data frame\cr
#'  cov.dat: \tab covariate data frame [or NA]
#'  }
#' @note Data files should be in PLINK ped/map format
#' @note Covariates file should have one column for subject
#' @note Do not include ped/map extension
#' @note Do include covariate extension IDs followed by one column per covariate
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @examples
#' init('/path/to','data')
#' init('/path/to','data','/path/to','covariates.csv')
init <- function(data.path,data.fname,cov.path=NA,cov.fname=NA) {
  # Load data from PLINK files and initialize the rare variant data frame

  # Load map, raw, and covariate data
  data <- load.data(data.path,data.fname,cov.path,cov.fname)
  map.dat <- data[[1]]
  raw.dat <- data[[2]]
  cov.dat <- data[[3]]
  
  # Identify rare variants
  # rs.id | chr | pos | minor.allele | allele.frequency
  rv.dat <- rare.vars(data.path,data.fname,map.dat)
  
  return(list("rv.dat"=rv.dat,"raw.dat"=raw.dat,"cov.dat"=cov.dat))
}

## ---------------------------------------------- ##
##               Support Functions                ##
## ---------------------------------------------- ##

rare.vars <- function(data.path,data.fname,map.dat) {
  # Identify the major/minor allele in all SNPs
  
  snp.dat <- allele.frequency(data.path,data.fname)
  snp.dat$SNP <- sapply(snp.dat$SNP,function(x){sub("[[:punct:]]",".",x)})

  # Separate the common and rare variants
  rv.dat <- rare.snps(snp.dat)
  
  # Use the map file to add position information
  rv.dat <- add.pos(rv.dat,map.dat)

  return(rv.dat)
}

allele.frequency <- function(data.path,data.fname) {
  # Given a PED/MAP path and basename, use PLINK to
  # return the rsIDs with frequency > 0 and < 0.05
  
  # Build the full file
  fin  <- paste(data.path,'/',data.fname,sep='')
  
  # Call PLINK to generate frequency file
  system(paste("plink --noweb --file",fin,"--freq --allow-no-sex",sep=' '))
  
  # Give PLINK time to write the file
  Sys.sleep(2)
  
  # Read the frequency file
  freq.dat <- read.table("plink.frq",sep='',header=TRUE,col.names=c("CHR","SNP","MA","MajA","MAF","NCHROBS"))
    
  # Subset the data to include only SNP, Minor Allele, and Frequency info
  freq.dat <- subset(freq.dat,select=c("SNP","MA","MAF"))
}

rare.snps <- function(freq.dat) {
  # Removes all common variants and returns the rare subset
  rare.snp.dat <- subset(freq.dat,MAF>0 & MAF<0.05)
}

add.pos <- function(rv.dat,map.dat) {
  # Uses the MAP file to add position info to the rare variant data frame
  
  col.names <- append(names(rv.dat),c("CHR","POS"))
  rv.dat <- merge(rv.dat,map.dat,by="SNP")
  rv.dat <- subset(rv.dat,select=col.names)
}

load.data <- function(data.path,data.fname,cov.path,cov.fname) {
  
  # Initialize all file names
  ped.file  <- paste(data.path,"/",data.fname,".ped",sep='')
  map.file  <- paste(data.path,"/",data.fname,".map",sep='')
  raw.file  <- paste(data.path,"/",data.fname,".raw",sep='')
  cov.file = NA
  if (!is.na(cov.path) & !is.na(cov.fname)) {
    cov.file <- paste(cov.path,"/",cov.fname,".cov",sep='')
  }
  
  # If the RAW file is missing, use PLINK to create it
  if (!file.exists(raw.file)) {
    system(paste("plink --noweb --file",paste(data.path,"/",data.fname,sep=''),"--recodeA --out",paste(data.path,"/",data.fname,sep='')))
  }
  
  # Load the data frames
  map.dat  <- read.table(map.file,sep=' ',header=FALSE,col.names=c('CHR','SNP','DIST','POS'))
  map.dat$SNP <- sapply(map.dat$SNP,function(x){sub("[[:punct:]]",".",x)})
  raw.dat  <- read.table(raw.file,sep=' ',header=TRUE)
  names(raw.dat) <- sapply(names(raw.dat),function(x){x <- strsplit(x,'_')[[1]][1]; sub("[[punct:]]",".",x)})
  cov.dat = NA
  if (!is.na(cov.file)) {
    cov.dat <- read.table(cov.file,sep='',header=TRUE)
  }
  
  return(list(map.dat,raw.dat,cov.dat)) 
}