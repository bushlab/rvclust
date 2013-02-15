# Project: RVCLUST
# File     clustering.r
#
# Author: R. Michael Sivley
#         Center for Human Genetics Research
#         Vanderbilt University Medical Center
# Email:  mike.sivley@vanderbilt.edu
#
# Description:
#
## ------------------------------------------------------------------------- ##

create.cluster.dat <- function(k,rv.dat,states.dat=NA,fitness=NA) {
  # Uses the clusters array, rare variant data frame, and chomratin states
  # to generate a new cluster.data frame
  
  # Initialize with the cluster IDs
  cluster.dat <- data.frame(CLUSTERID=1:k)
  
  # Add the sum of squares for each cluster
  if (!is.na(fitness)) {
    cluster.dat$FIT <- fitness
  }
  
  # Add the state information
  cluster.dat$STATE <- sapply(cluster.dat$CLUSTERID,function(x,dat){dat[dat$CLUSTERID==x,]$STATE[1]},dat=rv.dat)
  # Relabel the states with their groups
  label.states <- function(x) {if (x==1) {'promoter'} else if (x==2) {'enhancer'} else if (x==3) {'transcription'} else if (x==4) {'insulator'} else {'insulator'}}
  cluster.dat$STATE <- sapply(cluster.dat$STATE,label.states)
  
  # Add the chromosome information
  cluster.dat$CHR <- sapply(cluster.dat$CLUSTERID,get.chr,rv.dat=rv.dat)
  
  # Determine the minimum base position
  cluster.dat$MIN.BP <- sapply(cluster.dat$CLUSTERID,get.min.bp,rv.dat=rv.dat)
  
  # Determine the maximum base position
  cluster.dat$MAX.BP <- sapply(cluster.dat$CLUSTERID,get.max.bp,rv.dat=rv.dat)

  # Determine the size of each cluster
  
  cluster.dat$SIZE   <- sapply(1:k,function(i){nrow(rv.dat[rv.dat$CLUSTERID==i,])})
  
  return(cluster.dat)
  
  # Determine the chromatin state distribution
  #dists <- sapply(cluster.dat$CLUSTERID,get.chrm.dist,rv.dat=rv.dat,states.dat=states.dat)
}

collapse.clusters <- function(rv.dat,raw.dat,burden,column.only=FALSE) {
  # Collapse the raw dataset by cluster
  
  if (class(rv.dat) == "list") {
    rv.dat <- rv.dat[[1]]}
  
  # Extract the unique cluster IDs
  clusterids <- as.character(unique(rv.dat$CLUSTERID))
  clusterids <- sort(clusterids)
  
  # Use the clusterids and rv.dat to partition raw.dat
  # into clusters and collapse
  collapsed.vec <- sapply(clusterids,collapse.cluster,rv.dat=rv.dat,raw.dat=raw.dat,burden=burden)
  
  # If column.only is specified, then only one cluster was passed
  # and only that collapsed column should be returned
  if (column.only) {
    return(collapsed.vec[,1])}
  
  # Convert the collapsed data vector into a cluster-wise data frame
  else {
    #red.raw.dat <- subset(raw.dat,select=c('FID','IID','PAT','MAT','SEX','PHENOTYPE'))
    red.raw.dat <- raw.dat[,names(raw.dat) %in% c('FID','IID','PAT','MAT','SEX','PHENOTYPE')]
    collapsed.dat <- data.frame(collapsed.vec)
    names(collapsed.dat) <- clusterids
    collapsed.dat <- append(red.raw.dat,collapsed.dat)
    return(collapsed.dat)
  }
}

collapse.cluster <- function(id,rv.dat,raw.dat,burden) {
  # Collapse a single cluster of SNPs according to the burden flag
  # If burden is true, cluster value is the sum of minor alleles across member SNPs
  # If burden is false, cluster value is 1 if any minor allele present, 0 otherwise
  
  allClass <- function(x) {unlist(lapply(unclass(x),class))}
  
  # Reduce to SNPs in this cluster (id)
  #rv.clustered.dat <- subset(rv.dat,CLUSTERID==id)
  rv.clustered.dat <- rv.dat[rv.dat$CLUSTERID==id,]
  
  # Extract the SNP names
  snps <- as.character(rv.clustered.dat$SNP)
  
  # Extract those SNP columns from raw.dat
  #snp.cols <- subset(raw.dat,select=names(raw.dat) %in% snps)
  #snp.cols <- raw.dat[snps]
  snp.cols <- raw.dat[,names(raw.dat) %in% snps]
  
  all.classes <- allClass(snp.cols)
  
  # Compute the row-wise sum of those columns as a vector
  val <- NA
  if (class(snp.cols) != "data.frame") {
    val <- snp.cols}
  else {
    val <- rowSums(snp.cols,na.rm=TRUE)}
  
  # If burden testing is not specified, convert to 0/1 coding
  if (!burden) {
    val <- (val > 0)*1
  }
  return(val)
}

get.min.bp <- function(id,rv.dat) {
  # Given a CLUSTERID and rare variant data frame, determine
  # the minimum base position in the cluster
  min.bp <- min(subset(rv.dat,CLUSTERID==id)$POS)
}

get.max.bp <- function(id,rv.dat) {
  # Given a CLUSTERID and rare variant data frame, determine
  # the maximum base position in the cluster
  max.bp <- max(subset(rv.dat,CLUSTERID==id)$POS)
}

get.chr <- function(id,rv.dat) {
  # Given a CLUSTERID and rare variant data frame, determine
  # the cluster's chromosome
  chr <- subset(rv.dat,CLUSTERID==id)$CHR[1]
}