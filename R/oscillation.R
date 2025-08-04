
## DO WAVELET ANALYSIS OF CELLS SORTED BY PHASE

#' Wavelet analysis of counts along a circular pseudophase.
#' @param phi a circular pseudophase assignement for the counts vector
#'     or matrix.
#' @param counts a vector or matrix of counts of the same dimension
#'     and order as the pseudophase phi.
#' @param ID if counts is a matrix, ID specifies the row number or row
#'     name for which the analysis is to be performed.
#' @param verb verbosity level, set to >2 to also get verbose output
#'     from the wavelet function.
#' @param ... arguments to \code{\link[WavletComp]{analyze.wavelet}}.
#'@export
get_wavelet <- function(phi, counts, ID, verb=1, ...) {

    ## default: expect a single counts vector,
    ## if ID is provided get data from full count matrix
    
    if ( !missing(ID) ) counts <- counts[ID,]

    ## ORDER
    if ( any(diff(phi)<0) ) {
        counts <- counts[order(phi)]
        phi <- sort(phi)
        if ( verb>0 )
            cat(paste("NOTE: re-ordered data along phi\n"))
    }

    ## time in fraction of cycle in radian
    dt <- mean(diff(phi)) # 
    
    ## pad circular data to remove edge effects
    
    n <- length(counts)
    n2 <- floor(n/2)
    counts <- c(counts[(n2+1):n],counts,counts[1:n2])
    ctme <- c(phi[(n2+1):n]-2*pi, phi, phi[1:n2]+2*pi)
    
    ## remove rare duplicated time points
    dups <- which(duplicated(ctme))
    if ( length(dups)>0 ) {
        warning(paste('removing', length(dups),
                      'duplicated time points\n'))
        counts <- counts[-dups]
        ctme <- ctme[-dups]
    }
    
    ## maximal period: take all
    upperPeriod <- length(counts)*dt
    
    clet <- WaveletComp::analyze.wavelet(data.frame(x=counts), "x",
                                         lowerPeriod = 2*dt, 
                                         upperPeriod = upperPeriod, 
                                         dt = dt, 
                                         ##loess.span = LOESS.SPAN,
                                         ##make.pval = NPERM>0,
                                         ##method = "white.noise",
                                         ##n.sim = NPERM, 
                                         verbose = verb>1,
                                         ...)
    ## RESULT: expanded/filtered  data 
    list(wlet=clet, data=cbind(ctme, counts))
}

## Wavelet coherence analysis between two counts vectors along a
## circular pseudophase.
## @param phi a circular pseudophase assignment for the counts matrix.
## @param counts a vector or matrix of counts where
##     and order as the pseudophase phi.
## @param ID1 if counts is a matrix, ID specifies the row number or
##     row name for which the analysis is to be performed.
## @param ID2 if counts is a matrix, ID specifies the row number or
##     row name for which the analysis is to be performed.
## @param verb verbosity level, set to >2 to also get verbose output
##     from the wavelet function.
## @param ... arguments to \code{\link[WavletComp]{analyze.wavelet}}.
##@export - TODO: implement for ID vector
get_coherence <- function(phi, counts, ID1, ID2, verb=0, ...) {


    clet <- analyze.coherency(df,
                              my.pair = c(iID, jID),
                              lowerPeriod = 2*dt, 
                              upperPeriod = upperPeriod, 
                              dt = dt, dj = 1/20,
                              loess.span = LOESS.SPAN,
                              make.pval = CNPERM>0,
                              method = "white.noise",
                              n.sim = CNPERM,
                              verb = verb>1)

}
