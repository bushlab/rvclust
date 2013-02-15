# Project: RVCLUST
# File     entropyClust.r
#
# Author: R. Michael Sivley
#         Center for Human Genetics Research
#         Vanderbilt University Medical Center
# Email:  mike.sivley@vanderbilt.edu
#
# Description:
#
## ------------------------------------------------------------------------- ##

# This file should be used to compute the optimal partitions
# for rare variants using position and chromatin state as the
# primary clustering variables and then optimizing to increase
# statistical power by maximizing entropy.

# Assuming the existence of a collapse.cluster function
# collapse.clusters <- function(vars.dat,raw.dat,burden)
#source('rvcluster.R')

# vec is vector of 0/1 values
# returns the frequency of 1
entropy <- function(vec) {
  
  # Compute the entropy of this binary vector
  freq <- length(vec[vec==1]) / length(vec)

  # If there were no positive cases, return an entropy of -1
  if (freq < 0.01) {return(-1)}

  # Otherwise, find the difference from 0.5, and subtract from 1
  # for the fitness [0.5..1]
  return(1-abs(0.5 - freq))
}

# vars.dat is a data.frame with only the desired columns
# vars.dat will be clustered over these columns (pos,state)
# data will be used to collapse the columns and compute entropy
entropyClust <- function(vars.dat,label.col,data,objfun,depth=0) {

  # If this is a recursive call and the data frame is embedded
  # in a list, then extract before processing
  if (class(vars.dat) == "list") {
    vars.dat <- vars.dat[[1]]}

  # Determine the objective function
  #FIXME: not sure why this isn't working, but it isn't a priority
  #FUN <- match.fun(objfun)
  FUN <- entropy

  # Exclude certain columns during clustering
  exclude <- c(label.col,"CLUSTERID")

  # If there are two or fewer rows left in vars.dat
  # (two columns in data), return as a leaf cluster
  if (nrow(vars.dat) <= 2) {
    #return(c(vars.dat))}
    return(list(vars.dat$SNP,NA,NA,-2))}

  # Collapse the collapsed parent for comparison
  vars.dat$CLUSTERID <- rep(1,nrow(vars.dat))
  parent.collapsed <- collapse.clusters(vars.dat,data,FALSE,column.only=TRUE)
  parent.fitness <- FUN(parent.collapsed)
  
  # If the parent's fitness is -1 (no positive cases), return as leaf
  if (parent.fitness == -1) {
    #return(c(vars.dat))}
    return(list(vars.dat$SNP,NA,NA,parent.fitness))}
 
  # Determine the optimal clustering of this partition
  # according to proximity and chromatin state
  vars.dat$CLUSTERID <- cluster::pam(vars.dat[,!(names(vars.dat) %in% exclude)],k=2,cluster.only=TRUE)
  vars.df.matrix <- split(vars.dat,vars.dat$CLUSTERID)
 
  # Determine the fitness of the children
  children.collapsed <- lapply(vars.df.matrix,collapse.clusters,raw.dat=data,burden=FALSE,column.only=TRUE)
  fitness.vec <- sapply(children.collapsed,FUN)
  
  # If both children have a fitness of -1 (no positive cases), return the parent as a leaf
  if (length(fitness.vec[fitness.vec == -1]) == length(fitness.vec)) {
    #leaf.vec <- c(vars.dat)
    return(list(vars.dat$SNP,NA,NA,parent.fitness))}

  # If either of the children are better fit (higher score) than the parent, continue
  if (length(fitness.vec[fitness.vec >= parent.fitness]) > 0) {
    children <- lapply(vars.df.matrix,entropyClust,label.col=label.col,data=data,objfun=objfun,depth=depth+1)
    leaf <- list(vars.dat$SNP,children[[1]],children[[2]])}

  # If the parent is better fit than either child, return the parent as a leaf
  else {
    #print('Ending branch. Returning partition as leaf.')
    #leaf.vec <- c(vars.dat)}
    leaf <- list(vars.dat$SNP,NA,NA,parent.fitness)}
  
  return(leaf)
}
