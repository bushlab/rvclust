#' Constructor for rvclustobject
#'
#'   Constructor for rvclustobject. Loads ped/map data from the\cr
#' specified path and filename. If a covariate path and\cr
#' filename are specified, the covariates will be loaded as well.\cr
#' If annotations are specified, rvclust will use GWAR to include\cr
#' the requested annotations. Burden testing can be enabled by\cr
#' setting the burden flag to true. A threshold can also be set\cr
#' such that only clusters with fitness exceeding the threshold\cr
#' are tested for association.
#'   The rvclust package contains wrappers for various clustering\cr
#' methods and statistical tests. Documentation on these wrappers\cr
#' is available. Simply pass the rvclustobject to any method with\cr
#' an rvclust wrapper and an updated rvclustobject will be returned.
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @method
#' @param pedmap.path path to pedmap files
#' @param pedmap.fname basename for pedmap files (no extension)
#' @param cov.path path to covariates file [optional]
#' @param cov.fname filename for covariates file (with extension) [optional]
#' @param annotations list of desired annotations [default ALL]
#' @param burden boolean indicating whether to burden test [default FALSE]
#' @param min.fit minimum cluster fitness to test (0..1) [default 0.0]
#' @return rvclustobject containing\cr
#'  \tabular{ll}{
#'  data: \tab ped, map, and covariate data
#'  variants: \tab rare variant data.frame\cr
#'  clusters: \tab empty list of clusters
#'  clusterinfo: \tab empty list of cluster info
#'  }
#' @note Data must be in PEDMAP format
#' @note Covariates - One column for subject; One column per covariate
#' @examples
#' rvclustobject('/path/to','pedmap')
#' rvclustobject(...,cov.path='/path/to',cov.fname='covariates.csv')
#' rvclustobject(...,annotations=c("CHROMATIN","OTHER"))
#' rvclustobject(...,burden=TRUE)
#' rvclustobject(...,min.fit=0.8)
rvclustobject <- function(pedmap.path,pedmap.fname,cov.path=NA,cov.fname=NA,
                annotations=NA,burden=FALSE,min.fit=0.0) {

	# Load map, raw, and covariate data
  data <- load.data(pedmap.path,pedmap.fname,cov.path,cov.fname)
  map.dat <- data[[1]]
  raw.dat <- data[[2]]
  cov.dat <- data[[3]]
  
  # Identify rare variants
  rv.dat <- rare.vars(pedmap.path,pedmap.fname,map.dat)
  
  # Create the rvclustobject and specify its class
  rv <- list("data"=list("ped"=raw.dat,"map"=map.dat,"cov"=cov.dat),
  		"variants"=rv.dat,"clusters"=NA,"clusterinfo"=NA,
      "annotations"=annotations,"burden"=burden,"min.fit"=min.fit)
  class(rv) <- "rvclustobject"

  return(rv)
}

#' S3 method definitions for rvclustobject
show <- function(x) {UseMethod("show",x)}
show.rvclustobject <- function(x) {print(c("data","variants","clusters","clusterinfo"))}

print <- function(x) {UseMethod("print",x)}
print.rvclustobject <- function(x) {print(c("data","variants","clusters","clusterinfo"))}

summary <- function(x) {UseMethod("summary",x)}
summary.rvclustobject <- function(x) {
  print("rvclust object:")
  print("Data:")
  print(summary(x$ped))
  print(summary(x$map))
  print(summary(x$cov))
  print("Variants:")
  print(summary(x$variants))
  print("Clusters:")
  print(summary(x$clusters))
  print("Cluster Info:")
  print(summary(x$clusterinfo))
}


## ---------------------------------------------- ##
##               Support Functions                ##
## ---------------------------------------------- ##

rare.vars <- function(pedmap.path,pedmap.fname,map.dat) {
  # Identify the major/minor allele in all SNPs
  
  snp.dat <- allele.frequency(pedmap.path,pedmap.fname)
  snp.dat$SNP <- sapply(snp.dat$SNP,function(x){sub("[[:punct:]]",".",x)})

  # Separate the common and rare variants
  rv.dat <- rare.snps(snp.dat)
  
  # Use the map file to add position information
  rv.dat <- add.pos(rv.dat,map.dat)

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

load.data <- function(pedmap.path,pedmap.fname,cov.path,cov.fname) {
  
  # Initialize all file names
  ped.file  <- paste(pedmap.path,"/",pedmap.fname,".ped",sep='')
  map.file  <- paste(pedmap.path,"/",pedmap.fname,".map",sep='')
  raw.file  <- paste(pedmap.path,"/",pedmap.fname,".raw",sep='')
  cov.file = NA
  if (!is.na(cov.path) & !is.na(cov.fname)) {
    cov.file <- paste(cov.path,"/",cov.fname,".cov",sep='')
  }
  
  # If the RAW file is missing, use PLINK to create it
  if (!file.exists(raw.file)) {
    system(paste("plink --noweb --file",paste(pedmap.path,"/",pedmap.fname,sep=''),"--recodeA --out",paste(pedmap.path,"/",pedmap.fname,sep='')))
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