
## GENERATE COHORT LISTS AND MATRIXES

#' Convert a gene cohort list or matrix to a vector of unique gene classes.
#' @export
cohort2clustering <- function(cohorts, n, na="na") {

    if ( inherits(cohorts, "list") ) {
        if ( missing(n) )  {
            warning("no total matrix size provided, using maximal index")
            n <- max(unlist(COHORTS))
        }
        cohorts <- get_cohorts(cohorts, n=nrow(genes))
    }
    if ( missing(n) ) n <- nrow(cohorts)
    
    cls <- rep(na, n)
    for ( j in 1:nrow(cohorts) ) {
        idx <- which(cohorts[j,]!=0)
        if ( any(cls[idx] != na) )
            stop("overlapping cohorts, can't be converted to clustering")
        cls[idx] <- rownames(cohorts)[j]
    }
    cls    
}

#' convert a gene clusterng vector to a gene cohort list
#' @export
clustering2cohorts <- function(cls, cls.srt) {

    if ( missing(cls.srt) ) cls.srt <- unique(cls)
    cohorts <- list()

    for ( cl in cls.srt )
        cohorts[[cl]] <- which(cls==cl)
    
    cohorts
}

## TODO:
## * provide this by mapping functions
## * optional cohorts: use row names to find cohorts in predefined lists,
#' Convert a gene cohort list to a boolean or size-normalized cohort matrix
#' @param cohorts a list of integer vectors, each providing the row indices
#' of the genes in the count table.
#' @param n the number of rows (genes) in the count table, if not provided
#' the maximal index will be used
#' @param normalize divide each boolean entry (1) by the cohort size
#' @export
get_cohorts <- function(cohorts, n, normalize=FALSE) {

    if ( missing(n) ) {
        warning("no total matrix size provided, using maximal index")
        n <- max(unlist(COHORTS))
    }
    
    ## generate boolean cohort matrix: cohort X genes
    cohm <- matrix(0, nrow=length(cohorts), ncol=n)
    rownames(cohm) <- names(cohorts)
    for ( j in seq_along(cohorts) ) 
        cohm[j, cohorts[[j]]] <- 1

    ## normalize by cohort size
    if ( normalize )
        cohm <- cohm/apply(cohm, 1, sum)

    cohm
}
