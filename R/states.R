
## CALCULATE THE COHORT STATE MATRIX

#' Normalize count table by total counts per cell
#' @param count table with genes in rows and cells in columns
#' @export
normalize_counts <- function(counts, check=TRUE) {

    total <- apply(counts, 2, sum, na.rm=TRUE)

    norm <- TRUE
    if ( check ) {
        ## check if it is already size-normalized
        if ( sd(total) < .Machine$double.eps ) {
            cat(paste("count table is already size-normalized: skipped\n"))
            return(counts)
        }
        if ( any(counts<0, na.rm=TRUE) ) {
            cat(paste("negative numbers detected, ",
                      "count table is normalized: skipped\n"))
            return(counts)
        }
    }
    ## TODO: faster way w/o double transpose
    counts <- t(t(counts)/total)

    ## adding total counts as attribute
    attr(counts, "total") <- total
    
    counts

}


## TODO: * automatic normalization to total gene counts * optional
## cohort recovery stats: how many are present?  * optional cohort
## optimization based on internal consistency

#' Generate a cohort expression state matrix.
#' @param counts a gene-wise read count table, with genes as rows and
#'     cells as columns.
#' @param cohorts cohort definition: a list where each list entry is a
#'     vector of indices in the count table.
#' @export
get_states <- function(counts, cohorts) {

    ## NOTE: cohorts is either already a normalized cohort matrix
    ## or a list of vectors with row indices of cohort genes in the count table
    
    if ( inherits(cohorts, "list" ) )
        cohorts <- get_cohorts(cohorts, n=nrow(counts), normalize=TRUE)
    colnames(cohorts) <- rownames(counts)

    ## handle NA counts
    counts[is.na(counts)] <- 0

    ## state =  cohorts X counts
    state <- cohorts %*% counts
    state
}
