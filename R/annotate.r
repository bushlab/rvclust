#' Applies annotations to variants for informed clustering
#'
#' ---------------------------------------------------------------------- '#\cr
#' Annotations allow RVCLUST to estimate which features are likely to have \cr
#' similar effects without looking at the outcomes in advance. This \cr
#' preserves much of the statistical integrity lost with supervised \cr
#' learning. These annotations may also aid in interpreting results. \cr
#' ---------------------------------------------------------------------- '#\cr
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @param rv rvclustobject
#' @return annotated rvclustobject
annotate <- function(rv,annotations=NA,file=NA) {
  rv.dat <- rv$variants
  rv$annotations <- annotations

  ## Chromatin State Annotation ##
  if ('CHROMATIN' %in% annotations | is.na(annotations)) {
    
    # NOTE: Package will make GWAR API calls when a byloc API
    # is available. Until that time, chromatin annotations are
    # loaded from a local file, and limited to cell line: 
    # wgEncodeBroadHmmGm12878HMM

    # Package chromatin flat file and read directly
    data(Histone_HMM_subset)
    chrom.table <- Histone_HMM_subset
    
    
    get.chrm.state <- function(x,states.dat) {
      # Given a base position, determine the chromatin state
      states <- states.dat[states.dat$chromosome==paste('chr',x[2],sep='') & states.dat$start <= x[1] & x[1] < states.dat$end,]$state[1]
      states <- sapply(states,compress.states)
    }
    
    # Annotate rv.dat with chromatin state
    states <- apply(cbind(rv.dat$POS,rv.dat$CHR),1,get.chrm.state,states.dat=chrom.table)
    rv.dat$CHROMATIN <- states
  }

  ## PDB Coordinate Annotation ##
  if ('PDB' %in% annotations) {
    
    pdbmap <- read.table(file,sep='',header=TRUE)

    get.pdb.annotations <- function(x,dim,pdbmap) {
      # Given a chromosome and base position, determine the chromatin state
      coord <- pdbmap[pdbmap$chr==x[2] & pdbmap$start <= x[1] & x[1] <= pdbmap$end,dim][1]
    }

    # Add all PDB annotations
    rv.dat$PDBID  <- apply(cbind(rv.dat$POS,rv.dat$CHR),1,get.pdb.annotations,dim='pdbid',pdbmap=pdbmap)
    rv.dat$CHAIN  <- apply(cbind(rv.dat$POS,rv.dat$CHR),1,get.pdb.annotations,dim='chain',pdbmap=pdbmap)
    rv.dat$SEQRES <- apply(cbind(rv.dat$POS,rv.dat$CHR),1,get.pdb.annotations,dim='seqres',pdbmap=pdbmap)
    rv.dat$PDB_x  <- apply(cbind(rv.dat$POS,rv.dat$CHR),1,get.pdb.annotations,dim='x',pdbmap=pdbmap)
    rv.dat$PDB_y  <- apply(cbind(rv.dat$POS,rv.dat$CHR),1,get.pdb.annotations,dim='y',pdbmap=pdbmap)
    rv.dat$PDB_z  <- apply(cbind(rv.dat$POS,rv.dat$CHR),1,get.pdb.annotations,dim='z',pdbmap=pdbmap)
    
    # If PDB is the only annotation, remove all variants without PDB annotations
    if (length(annotations) == 1) {
      rv.dat <- rv.dat[complete.cases(rv.dat),]
    }
  }
  
  rv$variants <- rv.dat
  return(rv)
}

get.chrm.state <- function(x,states.dat) {
  # Given a base position, determine the chromatin state
  states <- states.dat[states.dat$chromosome==paste('chr',x[2],sep='') & states.dat$start <= x[1] & x[1] < states.dat$end,]$state[1]
  states <- sapply(states,compress.states)
}

compress.states <- function(state) {
  # Compress states
  
  if (is.na(state)) {
    state=NA
  }
  else if (state=='1_Active_Promoter' | state=='2_Weak_Promoter' | state=='3_Poised_Promoter' | state=='12_Repressed') {
    #state='promoter'
    state = 1
  }
  else if (state=='4_Strong_Enhancer' | state=='5_Strong_Enhancer' | state=='6_Weak_Enhancer' | state=='7_Weak_Enhancer') {
    #state='enhancer'
    state = 2
  }
  else if (state=='9_Txn_Transition' | state=='10_Txn_Elongation' | state=='11_Weak_Txn') {
    #state='transcription'
    state = 3
  }
  else if (state=='8_Insulator') {
    #state='insulator'
    state = 4
  }
  else if (state=='13_Heterochrom/lo' | state=='14_Repetitive/CNV' | state=='15_Repetitive/CNV') {
    #state='nonactive'
    state = 5
  }
  state
}
