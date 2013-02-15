#' Uses GWAR to annotate genetic variants
#'
#' @export
#' @param rv.dat rare variant data frame
#' @param annotations character list of desired annotations, all by default
#' @return annotated rare variant data frame
#' @seealso \code{\link{init}}
#' @author R Michael Sivley \email{mike.sivley@@vanderbilt.edu}
#' @examples
#' init(rv.dat)
annotate <- function(rv,annotations=NA) {
  rv.dat <- rv$variants

  ## Chromatin State Annotation ##
  if ('chromatin' %in% annotations | is.na(annotations)) {
    
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
    rv.dat$STATE <- states
  }
  
  rv$variants <- rv.data
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