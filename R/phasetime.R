### PSEUDOPHASE IN TIME: ERGODIC RATE ANALYSIS

## Ergodic rate analysis (ERA, Kafri et al. 2013, Box 2): in an
## asynchronous population that grows exponentially at rate
## a=log(2)/T, the speed along the cycle is v = a*(2-F)/f, with f the
## density of cells and F its cumulative fraction from division.
## Integrating dt = dl/v gives the time since division, T*(1-log2(2-F)).
## Only the order of the cells enters, the density is in F. See
## ROADMAP.md for the derivation and assumptions.

#' Cell age from the fraction of cells since division (ERA).
#'
#' In an asynchronous population that grows exponentially with
#' doubling time \code{T}, a cell that has a fraction \code{F} of all
#' cells between division and itself has the age
#' \code{T*(1-log2(2-F))} (ergodic rate analysis, Kafri et al. 2013).
#' @param F fraction of cells between division and the cell, in 0:1.
#' @param T doubling time of the population; the default 1 gives the
#'     age as a fraction of the cycle.
#' @return the age, in units of \code{T}.
#' @seealso \link{era_fraction}, the inverse, and \link{phase_time}.
#' @export
era_age <- function(F, T=1) T*(1 - log2(2 - F))

#' Fraction of cells younger than an age (inverse of ERA).
#' @param age age since division, in units of \code{T}.
#' @inheritParams era_age
#' @return the fraction of cells \code{F} between division and
#'     \code{age}.
#' @export
era_fraction <- function(age, T=1) 2 - 2^(1 - age/T)

#' Time spent in a phase that starts at division, by group (ERA).
#'
#' A group of cells that receives a share \code{share} of the newborn
#' cells (both cells of each division count) and holds a fraction
#' \code{f} of all cells before it leaves the phase, stays in the phase
#' for \code{-T*log2(1-f/(2*share))}. With \code{share=1} (all cells)
#' this is \link{era_age}; in budding yeast, mothers and daughters
#' each get \code{share=1/2}.
#' @param f fraction of all cells that are in the group and in the
#'     phase.
#' @param share the share of the newborn cells that enter the group.
#' @inheritParams era_age
#' @return the duration, in units of \code{T}; \code{NaN} if
#'     \code{f>=2*share}.
#' @export
era_duration <- function(f, T=1, share=1) -T*log2(1 - f/(2*share))

## distance along the cycle from an origin, in 0:2*pi
phase_since <- function(phi, origin=0) (phi - origin) %% (2*pi)

## fraction of cells (phase dcell, from the origin) before positions d
phase_fraction <- function(dcell, d)
    findInterval(d, sort(dcell), left.open=TRUE)/length(dcell)

## the phase origin: a number, or the moving average peak of a cohort
get_origin <- function(phases, origin, phase, states, window) {

    if ( is.numeric(origin) ) return(unname(origin[1]))

    if ( !missing(states) ) {
        pk <- state_phase(states, phases$rotation.phase[[phase]],
                          window=window)
        return(unname(pk[origin]))
    }
    if ( !phase %in% colnames(phases$x.phase) )
        stop("cohort origin without states requires x.phase$", phase)
    phases$x.phase[origin, phase]
}

#' Add the cell age to a phases object (ERA).
#'
#' Counts the cells from the origin (division) in the direction of
#' increasing phase and converts each cell's fraction \code{F} into its
#' age with \link{era_age}. Run it after \link{revert} and
#' \link{shift}, which do not update the age.
#' @inheritParams revert
#' @inheritParams era_age
#' @param origin the phase of division: a number, or the name of a
#'     cohort, whose moving average peak along the phase is used
#'     (\link{state_phase}, requires \code{states}) or else its phase in
#'     \code{phases$x.phase}.
#' @param phase the cell phase column of \code{phases$rotation.phase}
#'     that orders the cells, e.g. \code{"phi"}, \code{"theta"} (same
#'     order) or \code{"lambda"} (\link{principal_phase}).
#' @param states the cohort state matrix used for the phases, required
#'     for a cohort origin at its peak.
#' @param window fraction of cells for the moving average in
#'     \link{state_phase}.
#' @param name name of the new age column; the fraction is added as
#'     \code{<name>_F}.
#' @return the phases object with the cell age and fraction in
#'     \code{rotation.phase}, the cohort age in \code{x.phase} (if it
#'     has the phase column) and the settings in \code{phases$era}.
#' @export
phase_time <- function(phases, T=1, origin=0, phase='phi', states,
                       window=.05, name='age') {

    origin <- get_origin(phases, origin, phase, states, window)

    ## cells: fraction from the origin, mid-rank
    d <- phase_since(phases$rotation.phase[[phase]], origin)
    F <- (rank(d) - .5)/length(d)
    phases$rotation.phase[[name]] <- era_age(F, T)
    phases$rotation.phase[[paste0(name, '_F')]] <- F

    ## cohorts
    if ( phase %in% colnames(phases$x.phase) ) {
        Fc <- phase_fraction(d, phase_since(phases$x.phase[[phase]], origin))
        phases$x.phase[[name]] <- era_age(Fc, T)
    }

    phases$era <- list(T=T, origin=origin, phase=phase, name=name)
    phases$processing <- c(phases$processing, paste0("era:", name))
    phases
}

#' Landmarks along the cycle: cohort peaks and points of maximal change.
#'
#' Cohort peaks are the maxima of the moving average of each cohort
#' state along the phase (\link{state_phase}). A point of maximal
#' change between two cohort peaks is the cell with the maximal
#' smoothed slope \code{d(theta-phi)/dphi} (\link{get_segments},
#' \code{method="inflection"}) on the arc from the first to the second
#' peak, where the cell angle changes fastest per cell.
#' @inheritParams phase_time
#' @param between a named list of cohort pairs, \code{c(from, to)}, for
#'     points of maximal change.
#' @param slope name of the \link{get_segments} \code{"inflection"}
#'     result to use; computed with \code{spar} if missing.
#' @param spar smoothing parameter for \link{get_segments}.
#' @return named phases of the landmarks, cohort peaks first.
#' @export
phase_landmarks <- function(phases, states, between, phase='phi',
                            window=.05, slope='slope', spar=.001) {

    ph <- phases$rotation.phase[[phase]]
    peaks <- state_phase(states, ph, window=window)

    if ( missing(between) ) return(peaks)

    dcol <- paste0(slope, '_dtheta')
    if ( !dcol %in% colnames(phases$rotation.phase) )
        phases <- get_segments(phases, method='inflection', names=slope,
                               spar=spar, plot=FALSE)
    dth <- phases$rotation.phase[[dcol]]

    mx <- sapply(between, function(x) {
        d <- phase_since(ph, peaks[x[1]])
        arc <- which(d <= phase_since(peaks[x[2]], peaks[x[1]]))
        if ( length(arc)==0 ) return(NA)
        ph[arc][which.max(dth[arc])]
    })
    c(peaks, mx)
}

#' Ages of landmarks and durations between them (ERA).
#' @inheritParams phase_time
#' @param landmarks named phases of landmarks, e.g. from
#'     \link{phase_landmarks}.
#' @param T,origin,phase doubling time, origin and phase column; by
#'     default those of the last \link{phase_time} call.
#' @return a data frame with one row per landmark, sorted by age: the
#'     landmark \code{phase}, the fraction of cells \code{F} between the
#'     origin and the landmark, its \code{age}, the next landmark
#'     (\code{to}, \code{"division"} after the last) and the
#'     \code{duration} until then, all in units of \code{T}.
#' @export
phase_durations <- function(phases, landmarks, T, origin, phase) {

    era <- phases$era
    if ( missing(T) ) T <- era$T
    if ( missing(origin) ) origin <- era$origin
    if ( missing(phase) ) phase <- era$phase
    if ( is.null(T) || is.null(origin) || is.null(phase) )
        stop("T, origin and phase are required without phase_time")

    d <- phase_since(phases$rotation.phase[[phase]], origin)
    F <- phase_fraction(d, phase_since(landmarks, origin))
    res <- data.frame(landmark=names(landmarks), phase=unname(landmarks),
                      F=F, age=era_age(F, T))
    res <- res[order(res$age),]
    res$to <- c(res$landmark[-1], "division")
    res$duration <- c(res$age[-1], T) - res$age
    rownames(res) <- NULL
    res
}

#' Time spent in a phase that starts at division, per group of cells.
#'
#' Splits the cells between the origin and the end of a phase by a
#' grouping (e.g. daughters and mothers in G1) and converts each
#' group's fraction of all cells into its duration with
#' \link{era_duration}.
#' @inheritParams phase_durations
#' @param group a factor or vector, one entry per cell; \code{NA} cells
#'     count for the total but not for any group.
#' @param to the phase at which the phase ends, e.g. a landmark.
#' @param share the share of newborn cells per group, named by group;
#'     default: equal shares.
#' @return a data frame with one row per group: the number of cells in
#'     the phase \code{n}, their fraction of all cells \code{f}, the
#'     \code{share} and the \code{duration} (units of \code{T}); a last
#'     row \code{"all"} for all cells (share 1).
#' @export
group_durations <- function(phases, group, to, share, T, origin, phase) {

    era <- phases$era
    if ( missing(T) ) T <- era$T
    if ( missing(origin) ) origin <- era$origin
    if ( missing(phase) ) phase <- era$phase

    d <- phase_since(phases$rotation.phase[[phase]], origin)
    inp <- d < phase_since(to, origin)

    group <- as.factor(group)
    if ( missing(share) )
        share <- setNames(rep(1/nlevels(group), nlevels(group)),
                          levels(group))

    n <- c(table(group[inp]), all=sum(inp))
    f <- n/length(d)
    s <- c(share[levels(group)], all=1)
    data.frame(group=names(n), n=as.vector(n), f=as.vector(f),
               share=as.vector(s), duration=era_duration(f, T, s))
}

#' Profile of a per-cell value along age.
#'
#' Bins cells along age (or any phase) and summarizes a value in each
#' bin, e.g. the log total count as a size proxy, whose spread should
#' narrow where size is controlled.
#' @param age the cell age (or phase).
#' @param y the value, one per cell.
#' @param breaks number of equal-width bins, or the bin breaks.
#' @param lim the range of the bins, for a number of \code{breaks}.
#' @return a data frame with one row per bin: its \code{mid}, the
#'     number of cells \code{n}, the \code{mean}, \code{sd},
#'     \code{median}, and the coefficient of variation \code{cv}.
#' @export
age_profile <- function(age, y, breaks=20, lim=range(age, na.rm=TRUE)) {

    if ( length(breaks)==1 )
        breaks <- seq(lim[1], lim[2], length.out=breaks+1)
    bin <- cut(age, breaks=breaks, include.lowest=TRUE)

    res <- do.call(rbind, lapply(split(y, bin), function(x) {
        x <- x[!is.na(x)]
        c(n=length(x), mean=mean(x), sd=sd(x), median=median(x))
    }))
    data.frame(mid=(breaks[-1] + breaks[-length(breaks)])/2, res,
               cv=res[,'sd']/res[,'mean'], row.names=NULL)
}

#' Split cells by a two-component Poisson mixture of cohort counts.
#'
#' Cohort read counts \code{x} of cells with total counts \code{total}
#' are modelled as Poisson with rate \code{r_c*total} in one of two
#' classes \code{c}, fitted by EM. Used to call, e.g., daughter cells
#' from the counts of daughter-specific genes.
#' @param x cohort read counts, one per cell (e.g. a state times the
#'     cohort size times the total).
#' @param total the total counts of each cell.
#' @param iter maximal number of EM iterations.
#' @param tol stop when the log likelihood changes less.
#' @return a list with the posterior probability of the high-rate class
#'     \code{p}, the class \code{high} (\code{p>.5}), the rates
#'     \code{rate} and the class weights \code{weight}.
#' @export
split_counts <- function(x, total, iter=200, tol=1e-8) {

    ## start at the lower and upper halves of the rate
    rt <- x/total
    hi <- rt > median(rt)
    r <- c(sum(x[!hi])/sum(total[!hi]), sum(x[hi])/sum(total[hi]))
    r[1] <- max(r[1], 1e-3*r[2])
    w <- c(.5, .5)
    ll0 <- -Inf
    for ( i in seq_len(iter) ) {

        ## E-step, on the log scale
        l1 <- log(w[1]) + stats::dpois(x, r[1]*total, log=TRUE)
        l2 <- log(w[2]) + stats::dpois(x, r[2]*total, log=TRUE)
        mx <- pmax(l1, l2)
        ll <- sum(mx + log(exp(l1-mx) + exp(l2-mx)))
        p <- 1/(1 + exp(l1 - l2))

        ## M-step
        w <- c(1-mean(p), mean(p))
        r <- c(sum((1-p)*x)/sum((1-p)*total), sum(p*x)/sum(p*total))

        if ( abs(ll - ll0) < tol ) break
        ll0 <- ll
    }
    list(p=p, high=p>.5, rate=r, weight=w, loglik=ll)
}
