#' Annotate generic
#' 
#' @export
annotate <- function(x) {UseMethod("annotate",x)}
#' Uses GWAR to annotate the genetic variants in an
#' rvclustobject.
#'
#' The Genome-Wide Annotation Repository (GWAR) is a\cr 
#' resource provided by the Bush Lab at Vanderbilt\cr
#' University's that offers genome-wide annotations.\cr
#' More information on GWAR can be found at:\cr
#' \url{http://gwar.mc.vanderbilt.edu/}\cr
#' \cr
#' This method uses GWAR to annotate the variants\cr
#' in an rvclustobject with the annotations requested in\cr
#' the annotations list.\cr
#' \cr
#' Annotations will be appended to the ped object\cr
#'
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @export
#' @method annotate rvclustobject
#' @param rv rvclustobject
#' @return annotated rvclustobject
#' @seealso \code{\link{rvclust}}
annotate.rvclustobject <- function(rv) {
  rv.dat <- rv$variants
  annotations <- rv$annotations

  ## Chromatin State Annotation ##
  if ('CHROMATIN' %in% annotations | is.na(annotations)) {
    
    con <- RMySQL::dbConnect(RMySQL::MySQL(), user="mike",password="thorntonwellslab",
                 dbname="ucsc", host="gwar-dev.mc.vanderbilt.edu")
    chrom.table <- RMySQL::dbReadTable(con,"Histone_HMM_subset")
    capture <- RMySQL::dbDisconnect(con)
    
    get.chrm.state <- function(x,states.dat) {
      # Given a base position, determine the chromatin state
      
      states <- states.dat[states.dat$chromosome==paste('chr',x[2],sep='') & states.dat$start <= x[1] & x[1] < states.dat$end,]$state[1]
      states <- sapply(states,compress.states)
    }
    
    # Annotate rv.dat with chromatin state
    states <- apply(cbind(rv.dat$POS,rv.dat$CHR),1,get.chrm.state,states.dat=chrom.table)
    rv.dat$CHROMATIN <- states
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