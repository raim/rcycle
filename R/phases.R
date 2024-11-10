
### COHORT STATE MATRIX AND PSEUDOPHASE CALCULATION 

## TODO:
## * cohort mapping function via gene and species IDs,
## * pseudophase validation routines:
##     - correct cohort order,
##     - external pseudophase,
##     - external phase classes.



## return the rank of an angle in radian

#' Returns the rank of a phase angle in circular coordinates.
#' @param align align the original phase angle and the rank phase
#'     angle at phase 0.
#' @param center center the aligned phase angle such that all phases are
#'     in -pi:pi, using \link{align_phase}.
#' @export
phase_rank <- function(theta, align=TRUE, center=TRUE) {

    phi <- rank(theta)
    
    ## norm rank between -pi and pi
    phi <- 2*pi* phi/max(phi) - pi

    ## align at phase 0
    if ( align )
        phi <- align_phase(phi, target=theta, center=center)

    phi
}

#' Align two phase vectors at 0.
#' @export
align_phase <- function(phi, target, center=TRUE) {

    dph <- approx(x=target, y=phi, xout=0)$y
    shift_phase(phi, dph, center=center)

}

## simple circular interpolation
## TODO: better use approxfun.circular?
approx_phase <- function(x, y, xout, ...) {

    if ( any(xout > max(x+2*pi) ) | any(xout < min(x-2*pi)) )
        stop('xout of range')
    
    x <- c(x-2*pi, x, x+2*pi)
    y <- rep(y, 3)
    approx(x=x, y=y, xout=xout, ...)
}


### PHASES OBJECT FUNCTIONS 

#' Revert all phases if order is wrong.
#'
#' This function compares the calculated order of states to the input
#' order of the state matrix. If the reverse order better matches the
#' input order, all phases are reverted.
#'
#' @param phases phases object as returned by \link{get_pseudophase}.
#' @export
revert <- function(phases,  verb=1) {

    rphase <- revert_phase(phi=phases$x$phi)
    
    if ( rphase ) {
        if ( verb>0 )
            cat(paste("\treverting all phases\n"))

        ids <- c('phi', 'theta')
        
        phases <- lapply(phases, function(x) {
            for ( id in ids ) {
                ## NOTE: grepping phase angles with suffix, e.g. theta.s
                idx <- grep(paste0('^',id), names(x), value=TRUE)
                for ( i in idx )
                    x[[i]] <- -x[[i]]
            }
            if ( all(c('phi','order') %in% names(x)) )
                x$order <- order(x$phi)
            x
        })

        phases$processing <- c(phases$processing, "reverted")
        
    } else if ( verb>0 )
        cat(paste("\tphases already in correct order\n"))
    
    phases
}

#' Shift all phases and angles in a phase object.
#' @inheritParams get_pseudophase
#' @inheritParams shift_phase
#' @inheritParams revert
#' @export
shift <- function(phases, dphi, align=FALSE, center=FALSE, verb=1) {
    

    if ( verb>0 )
        cat(paste("\tshifting phases by", dphi, "\n"))

    ## shift phi in all items, re-order and phase align theta if present

    ## TODO: instead shift main theta and re-calculate rank phase?
    
    phases <- lapply(phases, function(x) {
        if ( 'phi' %in% names(x) ) {

            ## shift rank phase
            x$phi <- shift_phase(x$phi, dphi, center=TRUE)

            ## re-order
            if ( 'order' %in% names(x) )
                x$order <- order(x$phi)

            ## align all theta at 0
            if ( align ){
                idx <- grep('^theta', names(x), value=TRUE)
                for ( id in idx ) 
                    x[[id]] <- align_phase(phi=x[[id]],
                                           target=x$phi, center=center)
            }
        }
        x
    })

    
    phases$processing <- c(phases$processing, paste0("shift:", dphi))

    phases 
}


#' Calculate the difference between input and calculated state order.
#'
#' The differences are calculated as Levenshtein distances, using a
#' simple circular extension of the base R \link{adist} function.
#' @inheritParams revert
#' @export
evaluate_order <- function(phases) {

    ## simply compare the order of cohort phases with the input order
    ## of the state matrix, reflect in row order of cohort phases in pca$x
    state_order_distance(reference=rownames(phases$x),
                         test=rownames(phases$x)[phases$x$order])
}

#' Calibrate phase to a period.
#' @param period a period of the transcription cycle, in units of
#'     time.
#' @param phase the phase ID in the phases object to use for calibration
#'     to the transcription cycle period.
#' @inheritParams revert
#' @export
calibrate <- function(phases, period, phase='phi') {

    ## convert phase to time
    phases <- lapply(phases, function(x) {
        if ( phase %in% names(x) ) 
            x$time <- x[[phase]]/(2*pi) * period
        x
    })

    phases
}

classify <- function(phases) {}

## TODO: * expancd prcomp class instead of defining an new class!  *
## better way to reverse, * validate directly by comparing to row
## order of state,

#' get the pseudophase from a cohort state matrix
#' @param states a cohort expression state table, with cohort mean
#'     counts as rows and cells as columns, as provided by
#'     \link{get_states}.
#' @param verb output verbosity level.
#' @export
get_pseudophase <- function(states,
                            center, revert.phase=FALSE, window=.05,
                            segments=FALSE, spar=1e-3,
                            classify=FALSE, validate=FALSE,
                            ##use.states=FALSE,
                            row.center=TRUE,  # required: rm option?
                            verb=1) {

    cstates <- states
    
    ## ROW-CENTERING cohorts
    if ( row.center )
        cstates <- cstates - apply(cstates,1,mean) 
        
    ## PCA of cells
    ## scale: bring to unit variance and do COLUMN-CENTERING
    ## equiv to eigen(cor(cstates))
    phases <- prcomp(cstates, scale.=TRUE) 
    

    ## CELL PSEUDOPHASE from loadings of PC1 vs. PC2
    X <- phases$rotation[,1]
    Y <- phases$rotation[,2]

    ## get cell phase angle
    theta <- atan2(Y, X)
    amp <- sqrt(Y^2 + X^2)

    ## phase: rank(theta) in radian
    ## aligned at 0
    phi <- phase_rank(theta, align=TRUE, center=TRUE)

    
    ## COHORT PSEUDOPHASE
    cX <- phases$x[,1]
    cY <- phases$x[,2]

    ## get cohort phase angle
    ctheta <- atan2(cY, cX)
    camp <- sqrt(cY^2 + cX^2)

    ## interpolate cohort phase to cell rank phase via angle
    cphi <- approx_phase(x=theta, y=phi, xout=ctheta)$y
    

    ## TODO: remove jumps in theta already here?


    ## FURTHER PROCESS PHASES
    
    ## phase shift to center at a selected cohort
    if ( !missing(center) ) {

        ## centering based on center cohort maximal expression
        dph <- state_phase(states=cstates, phase=theta, center=center,
                           window=window, verb=verb)
    
        ## centering based on cohort phase angle
        ## appears to work WORSE than centering based on cohort expression;
        ## TODO: test different centerings/sortings in comparison
        if ( FALSE ) 
            dph <- ctheta[which(rownames(cstates)==center)]
        
        if ( verb>0 )
            cat(paste("phase shift by", round(dph,3), "\n"))

        ## cell phase shift
        theta <- shift_phase(theta, dph, center=TRUE)

        ## re-calculate rank phase
        phi <- phase_rank(theta, align=TRUE, center=TRUE)

        ## cohort phase shift
        ctheta <- shift_phase(ctheta, dph, center=TRUE)
        cphi <- shift_phase(cphi, dph, center=TRUE)
    }

    
    ## reverse phases if necessary
    if ( revert.phase ) {
        
        ##rphase2 <- revert_phase_state(phi=phi, states=state,  window=window)
        rphase <- revert_phase(phi=cphi)#, states=state,  window=window)

        ##if ( rphase!=rphase2 ) warning('orderings differ', rphase, rphase2)
        
        if ( rphase ) {
            if ( verb>0 )
                cat(paste("\treverting phases\n"))
            theta <- -theta
            phi <- -phi

            ctheta <- -ctheta
            cphi <- -cphi
        }
    }

   
    ## validate if state order (rows) is reproduced by phase order
    ## TODO: remove or generalize use, using PCA or state order
    sord <- NULL
    if ( validate ) {

        sord <- get_state_order(states=states[, order(phi)],
                                window=window, names=TRUE)
        ldst <- state_order_distance(reference=rownames(states),
                                     test=sord)
        if ( verb>0 )
            cat(paste("Levensthein order distance:", ldst, "\n"))
    }

    ## cohort ordering
    ord <- order(phi)

 

    ## SEGMENTS
    ## TODO: separate params for inflection and segments,
    ## TODO: skip inflections, just use segments and expand, add slopes etc.
    ## TODO: use these as alignment points for phase shifts
    if ( segments ) {
        stop('not implemented')
    }


    ## EXTENT PCA CLASS

    ## add default summary for pca
    ## add % var explained
    ##eigenvalues <- pca$sdev^2
    ##phases$variance <- eigenvalues/sum(eigenvalues)
    phases$summary <- summary(phases)$importance

    ## add cell phases to cell eigenvalues
    phases$rotation <- cbind.data.frame(order=order(phi),
                                     phi=phi,
                                     theta=theta,
                                     amplitude=amp, 
                                     phases$rotation)
    ## add cohort phases to cohort ?values?
    phases$x <- cbind.data.frame(order=order(cphi),
                              phi=cphi,
                              theta=ctheta,
                              amp=camp,
                              phases$x)

    ## simple order vector and distance
    ## TODO: instead add order vectors to phases$x
    ##phases$order <- rownames(phases$x)[phases$x$order]
    ##phases$distance <- state_order_distance(reference=rownames(states),
    ##                                     test=phases$order)
 
    ## add simple classification via max of log2 ratio
    ## TODO: PCA-based classification?
    if ( classify ) {
        cls <- get_classes(states=states)
        phases$rotation <- cbind(phases$rotation, class=cls)
    }
      
    ## extent PCA class
    class(phases) <- append("phases", class(phases))
    phases
}



#' Classify cells by maximal state log2 fold change.
#' @export
get_classes <- function(states) {

    cls <- row.names(states)

    ## cohort log2 fold change over mean
    nrm <- log2(states/apply(states,1,mean))

    ## for each cell get cohort with maximal fold change
    ccls <- apply(nrm, 2, function(x) {
        idx <- which.max(x)[1] ## NOTE: just take first of multiple!
        if ( length(idx)==0) idx <- NA # may not be required due to [1] above
        idx})
    
    cls[ccls]
}

#' Establish whether the phases should be reverted wrt state order.
#' @param phi a vector of phase angles (in radian).
#' @param ... parameters passed to \link{revert_phase_state}.
#' @inheritParams get_pseudophase
#' @export
revert_phase <- function(phi, states, ...) {

    ## use state maxima-based approach
    if ( !missing(states) )
        return(revert_phase_state(states=states, phi=phi, ...))
    
    
    ord <- order(phi)
    return(sum(diff(rev(ord))<0) <= sum(diff(ord)<0))

}

## TODO: avoid states, only do on PC1/2 angles
#' Establish whether the phases should be reverted wrt state order.
#' @inheritParams get_pseudophase
#' @inheritParams revert_phase
#' @inheritParams get_state_order
#' @export
revert_phase_state <- function(states, phi, window=.05) {

    idx <- get_state_order(states=states[, order(phi)], window=window)


    ## return IF phases should be reverted
    return(sum(diff(rev(idx))<0) <= sum(diff(idx)<0))

}




#' calculate Levenshtein distance between state/test and a reference order
#' @export
state_order_distance <- function(states, test, reference,  ...) {

    ## NOTE: either test string or state matrix is required
    if ( missing(test) )
        test <-  get_state_order(states=states, ...)
    if ( missing(reference) )
        reference <- rownames(states)

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


#' Get the order of states (rows of the cohort state matrix).
#'
#' NOTE: cells (columns of the state matrix) must be ordered!
#' @param states an ORDERED cohort expression state table, with cohort
#'     mean counts as rows and cells as columns, as provided by
#'     \link{get_states}. The columns must be ordered for the moving
#'     average calculation.
#' @param window fraction of cells to be used for the moving average
#'     of states.
#' @param circular treat the data as circular in the moving average
#'     calculation.
#' @export
get_state_order <- function(states, window=.05,
                            sides=2, circular=TRUE, names=FALSE) {

    ## window size for moving average, default: 5% of cells
    n <- ncol(states)*window
    
    ## find maxima of a moving average for each cohort (state row)
    nord <- rep(NA, nrow(states))
    for ( k in 1:nrow(states) ) {
        
        mavg <- stats::filter(states[k, ], rep(1/n, n),
                              sides=sides, circular=circular)
        nord[k] <- which.max(mavg)

    }
    
    ## return order of state maxima
    idx <- order(nord)

    ## return state names instead of row order
    if ( names )
        idx <- rownames(states)[idx]

    idx
}


## re-order:
## center and revert, assuming that rows in the state matrix reflect order,
## and provided a center
## TODO:
## * return all state phases, if center is not provided.
## * BETTER PHASE REVERSAL, fuse with general order validation.

#' Get the phase of a cohort state.
#' @export
state_phase <- function(states, phase, center, window=0.05, verb=0) {

    
    n <- ncol(states)*window # window size: default 5% of all cells
    ord <- order(phase) 
    
    ## get state to center 0: the peak of this state
    if ( is.character(center) )
        center <- which(rownames(states)==center)
    
    ## get peak cell of expression for each cohort
    ## of a moving average over phase-sorted cells
    nord <- rep(NA, nrow(states))
    for ( k in 1:nrow(states) ) {

        ## TODO: align with direct filter use and sides option
        ##       in get_state_order
        mavg <- ma(states[k, ord], n=n, circular=TRUE)
        nord[k] <- which.max(mavg)
        
    }
    
    ## get phase where cohort `center` has it's moving average peak
    ph0 <-  phase[ord][nord[center]]

    return(ph0)
}



#' Shift a phase vector by a certain phase.
#' @param phi phase angles to be shifted by dphi.
#' @param dphi phase angle by which the phase phi in the phases object
#'     will be shifted.
#' @param center re-center the shifted theta (align is true)
#'     between -pi:pi.
#' @export
shift_phase <- function(phi, dphi, center=TRUE) {

    phi <- phi - dphi
    if ( center )
        phi <- center_phase(phi)
    phi
    
}

#' Center phases between -pi and pi.
#' @param phi phase angle (radian).
#' @export
center_phase <- function(phi) {

    phi[phi< pi] <- phi[phi< pi] + 2*pi
    phi[phi>=pi] <- phi[phi>=pi] - 2*pi
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
remove_jumps <- function(phi, idx, center=TRUE, verb=1) {

 
    ## detect JUMPS in ANGLE and shift
    if ( missing(idx) )
        idx <- detect_jumps(phi)

    if ( length(idx)>1 ) warning("more than one jump")
    
    if ( length(idx)>0 ) {

        
        idx <- idx[1]
        
        ## shift phases before or after jump
        
        if ( idx > length(phi)/2 ) # append phases after jump
            phi[(idx+1):length(phi)] <- phi[(idx+1):length(phi)] + 2*pi
        else # prepend phases 
            phi[1:idx] <-  phi[1:idx] - 2*pi

        if ( verb>0 )
            cat(paste("shifting phase angles to remove the jump\n"))

        ## shift to -pi:pi
        if ( center & idx > length(phi)/2)
            phi <- phi + -pi - min(phi)
        else if ( center & idx <= length(phi)/2)
            phi <- phi + pi - max(phi)  
    }
    phi
}

