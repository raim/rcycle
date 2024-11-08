
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


#' Revert all phases if order is wrong.
#'
#' This function compares the calculated order of states to the input
#' order of the state matrix. If the reverse order better matches the
#' input order, all phases are reverted.
#'
#' @param phases
#' @param verb
#' @export
revert <- function(phases,  verb=1) {

    pca <- attr(phases, 'pca')

    
    rphase <- revert_phase(phi=pca$x$phi)
    
    if ( rphase ) {
        if ( verb>0 )
            cat(paste("\treverting all phases\n"))

        ids <- c('phi', 'theta')
        
        pca <- lapply(pca, function(x) {
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

        ## revert state order
        if ( 'order' %in% names(pca) )
            pca$order <- rev(pca$order)

        ## distance to reference
        if ( 'distance' %in% names(pca) )
            pca$distance <- state_order_distance(reference=rownames(pca$x),
                                                 test=pca$order)
        
        pca$processing <- c(pca$processing, "reverted")
        attr(phases, 'pca') <- pca
        
    } else if ( verb>0 )
        cat(paste("\tphases already in correct order\n"))
    
    phases    
}

#' Shift all phases and angles in a phase object.
#' @export
shift <- function(phases, dphi, align=FALSE, center=FALSE, verb=1) {
    
    pca <- attr(phases, 'pca')

    if ( verb>0 )
        cat(paste("\tshifting phases by", dphi, "\n"))

    ## shift phi in all items, re-order and phase align theta if present

    ## TODO: instead shift main theta and re-calculate rank phase?
    
    pca <- lapply(pca, function(x) {
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


    ## re-order and re-calculate distance
    if ( 'order' %in% names(pca) )
        pca$order <- rownames(pca$x)[pca$x$order]
    if ( 'distance' %in% names(pca) )
        pca$distance <- state_order_distance(reference=rownames(pca$x),
                                             test=pca$order)
    
    pca$processing <- c(pca$processing, paste0("shift:", dphi))
    attr(phases, 'pca') <- pca

    phases
}


#' Calculate the difference between input and calculated state order.
#'
#' The differences are calculated as Levenshtein distances, using a
#' simple circular extension of the base R \link{adist} function.
#' @param phases the phase object.
#' @export
evaluate_order <- function(phases) {

    phases <- attr(phases, 'pca')

    ## simply compare the order of cohort phases with the input order
    ## of the state matrix, reflect in row order of cohort phases in pca$x
    
    state_order_distance(reference=rownames(phases$x),
                         test=rownames(phases$x)[phases$x$order])
  
}

#' Calibrate phase to a period.
#' @export
calibrate <- function(phases, period, phase='phi') {

    pca <- attr(phases, 'pca')

    pca <- lapply(pca, function(x) {
        if ( phase %in% names(x) ) 
            x$time <- x[[phase]]/(2*pi) * period
        x
    })

    attr(phases, 'pca') <- pca
    phases
}

classify <- function(phases) {}

## TODO: * expancd prcomp class instead of defining an new class!  *
## better way to reverse, * validate directly by comparing to row
## order of state,

#' get the pseudophase from a cohort state matrix
#' @param state a cohort expression state table, with cohort mean
#'     counts as rows and cells as columns.
#' @export
get_pseudophase <- function(states,
                            center, revert.phase=FALSE, window=.05,
                            segments=FALSE, spar=1e-3,
                            classify=FALSE, validate=FALSE,
                            ##use.states=FALSE,
                            add.loadings=TRUE, # required for PCA plots
                            row.center=TRUE,  # required: rm option?
                            verb=1) {

    cstates <- states
    
    ## ROW-CENTERING cohorts
    if ( row.center )
        cstates <- cstates - apply(cstates,1,mean) 
        
    ## PCA of cells
    ## scale: bring to unit variance and do COLUMN-CENTERING
    ## equiv to eigen(cor(cstates))
    pca <- prcomp(cstates, scale.=TRUE) 
    

    ## CELL PSEUDOPHASE from loadings of PC1 vs. PC2
    X <- pca$rotation[,1]
    Y <- pca$rotation[,2]

    ## get cell phase angle
    theta <- atan2(Y, X)
    amp <- sqrt(Y^2 + X^2)

    ## phase: rank(theta) in radian
    ## aligned at 0
    phi <- phase_rank(theta, align=TRUE, center=TRUE)

    
    ## COHORT PSEUDOPHASE
    cX <- pca$x[,1]
    cY <- pca$x[,2]

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

    ## collect results
    res <- "hallo" #data.frame(order=ord, phase=phi, amplitude=amp, angle=theta)


    ## add simple classification via max of log2 ratio
    if ( classify ) {
        cls <- get_classes(states=states)
        res <- cbind(res, class=cls)
    }

    ## add all eigenvectors (loadings of PCA)
    if ( add.loadings ) {
        ##res <- cbind(res, pca$rotation)
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
    evals <- pca$sdev
    names(evals) <- paste0("PC",1:length(evals))
    attr(res, "eigenvalues") <- evals^2

    ## add cohorts
    ## get order phase via approx!
    ## TODO: ends, currently dirty via rule=2; how can this be done better?
    ## use circular approx: approxfun.circular
    ##cphi <- approx(x=theta, y=phi, xout=ctheta, rule=2)$y
    ##attr(res, "cohorts") <- data.frame(phase=cphi, amp=camp,
    ##                                   angle=ctheta, pca$x)

    ## add order as attribute
    if ( !is.null(sord) ) {
        ##attr(res, "order") <- sord    
        ##attr(res, "distance") <- ldst
    }


   ## TODO: extent PCA class instead of defining new class
    if ( TRUE ) {

        ## add default summary for pca
        pca$summary <- summary(pca)$importance
        
        pca$rotation <- cbind.data.frame(order=order(phi),
                                         phi=phi,
                                         theta=theta,
                                         amplitude=amp, 
                                         pca$rotation)
        pca$x <- cbind.data.frame(order=order(cphi),
                                  phi=cphi,
                                  theta=ctheta,
                                  amp=camp,
                                  pca$x)
        ## add % var explained
        ##eigenvalues <- pca$sdev^2
        ##pca$variance <- eigenvalues/sum(eigenvalues)

        
        pca$order <- rownames(pca$x)[pca$x$order]
        pca$distance <- state_order_distance(reference=rownames(states),
                                             test=pca$order)
 
        
        ## extent PCA class
        class(pca) <- append("phases", class(pca))

        ## TEMPORARY
        
        ## temporary: add as attribute until everything works
        attr(res, 'pca') <- pca

        ## temporary: use pca-based order if validate was false
        if ( is.null(sord) ) {

            ##attr(res, "order") <- pca$order
            ##attr(res, "distance") <- pca$distance
       }
        
    }

    ## DEFINE A CLASS 
    class(res) <- append("phases", class(res))
    
    res
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
#' @export
revert_phase_state <- function(states, phi, window=.05) {

    ## reverse order if most are negative
    ## NOTE: this relies on state matrix row order reflecting actual order!
    ##       it should be optional and in a function, taking alternative orders
    ##       as argument.
    ## TODO: better way to compare orders, currently fails for proline data,
    ## NOTE: the problem is related to validation by order, however,
    ##       to fully validate we need to account for circular data
    ##idx <- order(nord)

    idx <- get_state_order(states=states[, order(phi)], window=window)


    ## return IF phases should be reverted
    return(sum(diff(rev(idx))<0) <= sum(diff(idx)<0))

    ## TODO: do this in an external function on the phases object/class!
    ## reverse phases if the reverse order is equally good
    ## TODO: convert order to circular coordinates and use circular diff,
    ## and simply include first as last
    ##if ( sum(diff(rev(idx))<0) <= sum(diff(idx)<0) ) {
    ##    if ( verb>0 )
    ##        cat(paste("reversing phases\n"))
    ##    phi <- -phi
    ##}
    ##phi
}




#' calculate Levenshtein distance between state/test and a reference order
#' @export
state_order_distance <- function(state, test, reference,  ...) {

    ## NOTE: either test string or state matrix is required
    if ( missing(test) )
        test <-  get_state_order(states=state, ...)
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

## TODO: fix this, do not shift both phase and angle by the same
## delta, but shift one and then interpolate the others to 0 alignment

#' Shift all phases in a phase object.
#' @export
shift_phases <- function(phases, dph, center=TRUE, verb=0) {

    ## cell phases
    if ( verb>0 ) cat(paste("shift cell phases\n"))

    ## shift main phase
    phases$phase <- shift_phase(phases$phase, dph, center=center)

    ##dph2 <- approx(x=phases$phase, y=phases$angle, xout=0)$y
    ##phases$angle <- shift_phase(phases$angle, dph2, shift=TRUE)

    phases$angle <- shift_phase(phases$angle, dph, center=center)
    phases$order <- order(phases$phase)
       
    ## cohort phases
    if ( "cohorts"%in%names(attributes(phases)) ) {
        if ( verb>0 ) cat(paste("shift cohort phases\n"))
        coh <- attr(phases, "cohorts")
        coh$phase <- shift_phase(coh$phase, dph, center=center)
        coh$angle <- shift_phase(coh$angle, dph, center=center)
        coh$order <- order(coh$phase)
        attr(phases, "cohorts") <- coh
    }
    
    ## phases in segments and inflections
    if ( "segments"%in%names(attributes(phases)) ) {
        if ( verb>0 ) cat(paste("shift segment phases\n"))
        seg <- attr(phases, "segments")
        seg$x <- shift_phase(seg$x, dph, center=center)
        attr(phases, "segments") <- seg
    }
    if ( "inflections"%in%names(attributes(phases)) ) {
        if ( verb>0 ) cat(paste("shift inflection phases\n"))
        seg <- attr(phases, "inflections")
        seg$x <- shift_phase(seg$x, dph, center=center)
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

