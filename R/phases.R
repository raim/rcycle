
### COHORT STATE MATRIX AND PSEUDOPHASE CALCULATION 

## TODO:
## * cohort mapping function via gene and species IDs,
## * pseudophase validation routines:
##     - correct cohort order,
##     - external pseudophase,
##     - external phase classes.



## return the rank of an angle in radian ' returns the rank of a phase
#angle in circular coordinates
#' @param align align the original phase angle and the rank phase
#'     angle at phase 0.
#' @param shift shift the aligned phase angle such that all phases are
#'     in -pi:pi, using \link{phase_align}.
#' @export
phase_rank <- function(phi, align=TRUE, shift=TRUE) {

    rnk <- rank(phi)
    
    ## norm rank between -pi and pi
    rnk <- 2*pi* rnk/max(rnk) - pi

    ## align at phase 0
    if ( align ) rnk <- phase_align(rnk, target=phi, shift=shift)

    rnk
}

#' align to phase vectors at 0
#' @export
phase_align <- function(phi, target, shift=TRUE) {

    dph <- approx(x=target, y=phi, xout=0)$y
    phase_shift(phi, dph, shift=shift)

}

## TODO: * expancd prcomp class instead of defining an new class!  *
## better way to reverse, * validate directly by comparing to row
## order of state,
#' get the pseudophase from a cohort state matrix
#' @param state a cohort expression state table, with cohort mean
#'     counts as rows and cells as columns.
#' @export
get_pseudophase <- function(state,
                            center, revert.phase=FALSE, window=.05,
                            segments=FALSE, spar=1e-3,
                            classify=FALSE, validate=FALSE,
                            add.loadings=TRUE, # required for PCA plots
                            unit=FALSE, # convert PC1/2 to unit circle: remove
                            row.center=TRUE,  # required: rm option
                            verb=1) {

    cstate <- state
    
    ## ROW-CENTERING cohorts
    if ( row.center )
        cstate <- cstate - apply(cstate,1,mean) 
        
    ## PCA of cells
    ## scale: bring to unit variance and do COLUMN-CENTERING
    ## equiv to eigen(cor(cstate))
    pca <- prcomp(cstate, scale.=TRUE) 

    ## TODO: skip, this is better handled by order(phi)
    if ( unit ) 
        pca$rotation <- apply(pca$rotation, 2, minmax)*2 -1

    ## CELL PSEUDOPHASE from loadings of PC1 vs. PC2
    X <- pca$rotation[,1]
    Y <- pca$rotation[,2]

    ## get phase angle
    atn <- atan2(Y, X)
    amp <- sqrt(Y^2 + X^2)

    ## scaled phase: order in radian
    ## aligned at 0
    phase <- phase_rank(atn)

    
    ## COHORT PSEUDOPHASE
    cX <- pca$x[,1]
    cY <- pca$x[,2]
    catn <- atan2(cY, cX)
    camp <- sqrt(cY^2 + cX^2)


    ## FURTHER PROCESS PHASES:
    
    ## phase shift to center at a selected cohort
    if ( !missing(center) ) {

        ## centering based on center cohort maximal expression
        ph0 <- state_phase(state=cstate, phase=atn, center=center,
                           window=window, verb=verb)
    
        ## centering based on cohort phase angle
        ## appears to work WORSE than centering based on cohort expression;
        ## TODO: test different centerings/sortings in comparison
        if ( FALSE ) 
            ph0 <- catn[which(rownames(cstate)==center)]
        
        if ( verb>0 )
            cat(paste("phase shift by", round(ph0,3), "\n"))

        ## cell phase shift
        atn <- phase_shift(atn, ph0, shift=TRUE)

        ## re-calculate rank phase
        phase <- phase_rank(atn, align=TRUE, shift=TRUE)

        ## cohort phase shift
        catn <- phase_shift(catn, ph0, shift=TRUE)
    }

    
    ## reverse phases if necessary
    if ( revert.phase ) {
        rphase <- revert_phases(state, phase=phase,
                                window=window, verb=verb)
        if ( rphase ) {
            if ( verb>0 )
                cat(paste("\treverting phases\n"))
            atn <- -atn
            phase <- -phase

            catn <- -catn
        }
    }

    ## TODO: REMOVE JUMPS ALREADY HERE?
    ## and avoid downstream; problem:
    ##atn[order(phase)] <- remove_jumps(atn[order(phase)], shift=FALSE, verb=verb)
   
    ## validate if state order (rows) is reproduced by phase order
    sord <- NULL
    if ( validate ) {

        sord <- get_state_order(state[, order(phase)],
                                window=window, names=TRUE)
        ldst <- state_order_distance(reference=rownames(state),
                                     test=sord)
        if ( verb>0 )
            cat(paste("Levensthein order distance:", ldst, "\n"))
    }

    ## cohort ordering
    ord <- order(phase)

    ## collect results
    res <- data.frame(order=ord, phase=phase, amplitude=amp, angle=atn)


    ## add simple classification via max of log2 ratio
    if ( classify ) {
        cls <- get_classes(state=state)
        res <- cbind(res, class=cls)
    }

    ## add all eigenvectors (loadings of PCA)
    if ( add.loadings ) {
        res <- cbind(res, pca$rotation)
    }


    ## SEGMENTS
    ## TODO: separate params for inflection and segments,
    ## TODO: skip inflections, just use segments and expand, add slopes etc.
    ## TODO: use these as alignment points for phase shifts
    if ( segments ) {
        ## extrema of d(y-x)/dx
        segs <- get_segments(res, plot=FALSE,
                             ma.win=ceiling(nrow(res)*window),
                             spar=spar)
        attr(res, "segments") <- segs

        ## inflection points of y=f(x)
        infl <- get_inflections(res, plot=FALSE,
                                ma.win=ceiling(nrow(res)*window),
                                spar=spar)
        attr(res, "inflections") <- infl
    }

    ## add eigenvalues as attributes
    ## NOTE: "% variance explained" is eigenvalues/sum(eigenvalues)
    names(pca$sdev) <- paste0("PC",1:length(pca$sdev))
    attr(res, "eigenvalues") <- pca$sdev^2

    ## add cohorts
    ## get order phase via approx!
    ## TODO: ends, currently dirty via rule=2; how can this be done better?
    ## use circular approx: approxfun.circular
    cphase <- approx(x=atn, y=phase, xout=catn, rule=2)$y
    attr(res, "cohorts") <- data.frame(phase=cphase, amp=camp,
                                       angle=catn, pca$x)

    ## add order as attribute
    if ( !is.null(sord) ) {
        attr(res, "order") <- sord    
        attr(res, "distance") <- ldst
    }

    ## DEFINE A CLASS 
    class(res) <- append("phases", class(res))
    
    res
}

#' Classify cells by maximal state log2 fold change.
#' @export
get_classes <- function(state) {

    cls <- row.names(state)

    ## cohort log2 fold change over mean
    nrm <- log2(state/apply(state,1,mean))

    ## for each cell get cohort with maximal fold change
    ccls <- apply(nrm, 2, function(x) {
        idx <- which.max(x)[1] ## NOTE: just take first of multiple!
        if ( length(idx)==0) idx <- NA # may not be required due to [1] above
        idx})
    
    cls[ccls]
}


## TODO: avoid states, only do on PC1/2 angles
#' Establish whether the phases should be reverted wrt state order.
#' @export
revert_phases <- function(state, phase, window=.05, verb=1) {

    ## reverse order if most are negative
    ## NOTE: this relies on state matrix row order reflecting actual order!
    ##       it should be optional and in a function, taking alternative orders
    ##       as argument.
    ## TODO: better way to compare orders, currently fails for proline data,
    ## NOTE: the problem is related to validation by order, however,
    ##       to fully validate we need to account for circular data
    ##idx <- order(nord)

    idx <- get_state_order(state[, order(phase)], window=window)


    ## return IF phases should be reverted
    return(sum(diff(rev(idx))<0) <= sum(diff(idx)<0))

    ## TODO: do this in an external function on the phases object/class!
    ## reverse phases if the reverse order is equally good
    ## TODO: convert order to circular coordinates and use circular diff,
    ## and simply include first as last
    ##if ( sum(diff(rev(idx))<0) <= sum(diff(idx)<0) ) {
    ##    if ( verb>0 )
    ##        cat(paste("reversing phases\n"))
    ##    phase <- -phase
    ##}
    ##phase
}




#' calculate Levenshtein distance between state/test and a reference order
#' @export
state_order_distance <- function(state, test, reference,  ...) {

    ## NOTE: either test string or state matrix is required
    if ( missing(test) )
        test <-  get_state_order(state, ...)
    if ( missing(reference) )
        reference <- rownames(state)

    ## ABC-econding of classes in test and reference
    abc <- c(letters, toupper(letters))
    
    if ( length(unique(c(test,reference))) > length(abc) )
        stop("more classes than available letters")

    if ( !all(test %in% reference) )
        stop("not all test classes in reference class")
    
    ## convert to single letters
    ref <- abc[1:length(reference)]
    names(ref) <- reference

    tst <- ref[test]

    ref <- paste(ref,collapse="")
    tst <- paste(tst,collapse="")
    
    ## Levenshtein distance, simple circular
    dst <- Inf
    i <- 1
    n <- nchar(ref)
    tstst <- paste0(tst,tst)
    while( dst>0 & i<=n ) {
        stst <- substr(tstst, i, i + n - 1)
        i <- i+1
        ## calculate levenshtein distance and use if minimal
        dst <- min(c(dst, adist(ref, stst)[1,1]))
    }
    dst
}

## NOTE: CELLS (state columns) MUST BE ORDERED!
#' Get the order of states (rows of the cohort state matrix).
#' @export
get_state_order <- function(state, window=.05,
                            sides=2, circular=TRUE, names=FALSE) {

    ## window size for moving average, default: 5% of cells
    n <- ncol(state)*window
    
    ## find maxima of a moving average for each cohort (state row)
    nord <- rep(NA, nrow(state))
    for ( k in 1:nrow(state) ) {
        
        mavg <- stats::filter(state[k, ], rep(1/n, n),
                              sides=sides, circular=circular)
        nord[k] <- which.max(mavg)

    }
    
    ## return order of state maxima
    idx <- order(nord)

    ## return state names instead of row order
    if ( names )
        idx <- rownames(state)[idx]

    idx
}


get_phase <- function(phases) {


}

## re-order:
## center and revert, assuming that rows in the state matrix reflect order,
## and provided a center
## TODO:
## * return all state phases, if center is not provided.
## * BETTER PHASE REVERSAL, fuse with general order validation.

#' Get the phase of a cohort state.
#' @export
state_phase <- function(state, phase, center, window=0.05, verb=0) {

    
    n <- ncol(state)*window # window size: default 5% of all cells
    ord <- order(phase) 
    
    ## get state to center 0: the peak of this state
    if ( is.character(center) )
        center <- which(rownames(state)==center)
    
    ## get peak cell of expression for each cohort
    ## of a moving average over phase-sorted cells
    nord <- rep(NA, nrow(state))
    for ( k in 1:nrow(state) ) {
        mavg <- ma(state[k, ord], n=n, circular=TRUE)
        nord[k] <- which.max(mavg)
    }
    
    ## get phase where cohort `center` has it's moving average peak
    ph0 <-  phase[ord][nord[center]]

    return(ph0)
}

## TODO: fix this, do not shift both phase and angle by the same
## delta, but shift one and then interpolate the others to 0 alignment

#' Shift all phases in a phase object.
#' @export
shift_phases <- function(phases, dph, shift=TRUE, verb=0) {

    ## cell phases
    if ( verb>0 ) cat(paste("shift cell phases\n"))

    ## shift main phase
    phases$phase <- phase_shift(phases$phase, dph, shift=shift)

    ##dph2 <- approx(x=phases$phase, y=phases$angle, xout=0)$y
    ##phases$angle <- phase_shift(phases$angle, dph2, shift=TRUE)

    phases$angle <- phase_shift(phases$angle, dph, shift=shift)
    phases$order <- order(phases$phase)
       
    ## cohort phases
    if ( "cohorts"%in%names(attributes(phases)) ) {
        if ( verb>0 ) cat(paste("shift cohort phases\n"))
        coh <- attr(phases, "cohorts")
        coh$phase <- phase_shift(coh$phase, dph, shift=shift)
        coh$angle <- phase_shift(coh$angle, dph, shift=shift)
        coh$order <- order(coh$phase)
        attr(phases, "cohorts") <- coh
    }
    
    ## phases in segments and inflections
    if ( "segments"%in%names(attributes(phases)) ) {
        if ( verb>0 ) cat(paste("shift segment phases\n"))
        seg <- attr(phases, "segments")
        seg$x <- phase_shift(seg$x, dph, shift=shift)
        attr(phases, "segments") <- seg
    }
    if ( "inflections"%in%names(attributes(phases)) ) {
        if ( verb>0 ) cat(paste("shift inflection phases\n"))
        seg <- attr(phases, "inflections")
        seg$x <- phase_shift(seg$x, dph, shift=shift)
        attr(phases, "inflections") <- seg
    }
    
    if ( "time" %in% colnames(phases) ) {
        warning("recalculating time based on maximal range")
        phases$time <- (phases$phase)/(2*pi) * diff(range(phases$time))
    }
    
    ## return modified object
    phases
}

#' Shift a phase vector by a certain phase.
#' @export
phase_shift <- function(phi, dphi, shift=TRUE) {

    phi <- phi - dphi
    if ( shift ) {
        phi[phi< pi] <- phi[phi< pi] + 2*pi
        phi[phi>=pi] <- phi[phi>=pi] - 2*pi
    }
    phi
    
}


## TODO: more general solution?

#' Return the position of a jump
#' @export
detect_jumps <- function(py, max=-pi, verb=0) {
    which(diff(py) < max) 
}

#' Remove a jump in phases by shift all phases.
#' @param phi ordered list of phases
#' @export
remove_jumps <- function(phi, shift=TRUE, verb=1) {

 
    ## detect JUMPS in ANGLE and shift
    dph <- detect_jumps(phi)

    if ( length(dph)>1 ) warning("more than one jump")
    
    if ( length(dph)>0 ) {

        
        dph <- dph[1]
        
        ## shift phases before or after jump
        
        if ( dph > length(phi)/2 ) # append phases after jump
            phi[(dph+1):length(phi)] <- phi[(dph+1):length(phi)] + 2*pi
        else # prepend phases 
            phi[1:dph] <-  phi[1:dph] - 2*pi

        if ( verb>0 )
            cat(paste("shifting phase angles to remove the jump\n"))

        ## shift to -pi:pi
        if ( shift & dph > length(phi)/2)
            phi <- phi + -pi - min(phi)
        else if ( shift & dph <= length(phi)/2)
            phi <- phi + pi - max(phi)  
    }
    phi
}

