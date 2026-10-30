### DIVISION-CYCLE PHASE AND AGE: THE WHOLE CHAIN

## From counts to cell ages and phase durations: cohort states, the
## pseudophase from one set of cohorts (e.g. all periodic Cyclebase
## genes in peak-time bins), orientation and landmarks from another
## (e.g. the Spellman 1998 cohorts in cycle order), the ERA age, and
## noise flags. Plus confidence intervals of the durations by
## resampling cells, and time weights to turn snapshot (cell) averages
## into averages over the cycle.

#' Cohort states from counts.
#'
#' The state of a cohort in a cell is the mean over its genes of the
#' counts divided by the cell's total, as in \link{get_states} with
#' \link{normalize_counts}.
#' @param counts gene counts, genes x cells (also a sparse
#'     \code{Matrix}), with gene IDs as row names.
#' @param total total counts per cell; default: column sums.
#' @param cohorts named list of gene IDs per cohort; genes not in
#'     \code{counts} are ignored.
#' @return the state matrix, cohorts x cells.
#' @export
cohort_states <- function(counts, total, cohorts) {

    ## dim and products of sparse matrices need the Matrix methods
    if ( !is.matrix(counts) ) loadNamespace("Matrix")
    if ( missing(total) ) total <- Matrix::colSums(counts)
    w <- matrix(0, nrow=length(cohorts), ncol=nrow(counts),
                dimnames=list(names(cohorts), rownames(counts)))
    for ( k in names(cohorts) ) {
        g <- intersect(cohorts[[k]], rownames(counts))
        if ( length(g)==0 ) stop("no genes of cohort ", k, " in counts")
        w[k, g] <- 1/length(g)
    }
    st <- as.matrix(w %*% counts)
    st <- t(t(st)/total)
    colnames(st) <- colnames(counts)
    st
}

#' Division-cycle phase and cell age from counts.
#'
#' The whole chain: cohort states (\link{cohort_states}); the
#' pseudophase from \code{cohorts} (\link{get_pseudophase}); its
#' direction such that the peaks of the landmark cohorts \code{marks}
#' come in the given order; per-cell quality against Poisson noise
#' (\link{cell_quality}), optionally removing noise cells and
#' recomputing; the ERA age with the origin at the peak of the first
#' landmark cohort (\link{phase_time}); landmarks (cohort peaks and
#' points of maximal change between consecutive cohorts,
#' \link{phase_landmarks}) and their durations (\link{phase_durations}).
#' @inheritParams cohort_states
#' @param cohorts named list of gene IDs for the phase, in cycle order.
#' @param marks named list of gene IDs of the landmark cohorts, in cycle
#'     order, the first at division (e.g. M/G1).
#' @param T doubling time; default 1: ages in units of the cycle.
#' @param between pairs of landmark cohorts for points of maximal change;
#'     default: consecutive \code{marks} after the first, named by their
#'     concatenated names without prefix (e.g. \code{"G1S"}).
#' @param exclude gene IDs to leave out of the phase cohorts, e.g. the
#'     genes of a cohort that will be read against the phase.
#' @param B number of noise cells per cell (\link{cell_quality}); 0: no
#'     quality.
#' @param noise.rm remove noise cells and recompute the phase.
#' @param alpha noise threshold, see \link{cell_quality}.
#' @param spar smoothing for the points of maximal change
#'     (\link{get_segments}).
#' @param window fraction of cells for moving averages.
#' @param seed optional random seed.
#' @param verb verbosity.
#' @return a phases object (\link{get_pseudophase}) with the ERA age in
#'     \code{rotation.phase}, the landmark states \code{marks}, the
#'     \code{landmarks}, the \code{durations}, \code{order.ok} (landmark
#'     peaks in the given order), \code{between}, and \code{cells}, the
#'     columns of \code{counts} that were used.
#' @export
division_phase <- function(counts, total, cohorts, marks, T=1,
                           between, exclude, B=20, noise.rm=FALSE,
                           alpha=.05, spar=.001, window=.05, seed,
                           verb=0) {

    if ( missing(total) ) total <- Matrix::colSums(counts)
    if ( !missing(seed) ) set.seed(seed)
    if ( !missing(exclude) )
        cohorts <- lapply(cohorts, setdiff, exclude)
    m <- names(marks)
    if ( missing(between) ) {
        between <- lapply(2:(length(m)-1), function(i) m[i:(i+1)])
        names(between) <- sapply(between, function(x)
            paste(sub(".*_", "", x), collapse=""))
    }

    st <- cohort_states(counts, total, cohorts)
    mk <- cohort_states(counts, total, marks)
    n <- sapply(cohorts, function(g) length(intersect(g, rownames(counts))))

    P <- orient_phase(st, mk, window=window)
    cells <- seq_len(ncol(st))
    if ( B > 0 ) {
        P <- cell_quality(P, st, total, n, B=B, alpha=alpha)
        if ( noise.rm ) {
            cells <- which(!P$rotation.phase$noise)
            if ( verb>0 ) cat("removing", ncol(st)-length(cells),
                              "noise cells\n")
            q <- P$rotation.phase[cells, c("sd", "r2", "amplitude.raw",
                                           "noise.p", "noise.p.r2",
                                           "noise")]
            st <- st[, cells]
            mk <- mk[, cells]
            P <- orient_phase(st, mk, window=window)
            P$rotation.phase <- cbind(P$rotation.phase, q)
        }
    }

    P <- get_segments(P, method='inflection', names='slope', spar=spar,
                      plot=FALSE)
    P <- phase_time(P, T=T, origin=m[1], states=mk, window=window)
    lm <- phase_landmarks(P, mk, between=between, window=window)
    dur <- phase_durations(P, lm)

    P$marks <- mk
    P$between <- between
    P$landmarks <- lm
    P$durations <- dur
    P$order.ok <- identical(dur$landmark[dur$landmark %in% m], m)
    P$cells <- colnames(counts)[cells]
    P
}

## pseudophase, directed such that the landmark peaks come in order
orient_phase <- function(st, mk, window=.05) {
    P <- get_pseudophase(st, revert.phase=TRUE, verb=0)
    m <- rownames(mk)
    dst <- sapply(c(1, -1), function(s) {
        pk <- state_phase(mk, s*P$rotation.phase$phi, window=window)
        state_order_distance(reference=m,
                             test=names(sort(phase_since(pk, pk[m[1]]))))
    })
    if ( dst[2] < dst[1] ) P <- revert(P, force=TRUE, verb=0)
    P
}

#' Confidence intervals of landmark ages and durations.
#'
#' Resamples the cells with replacement and recomputes the origin
#' (peak of the first landmark cohort), the landmarks and their ERA
#' ages; the smoothed slope for the points of maximal change is kept
#' from the full data. Requires a result of \link{division_phase}.
#' @inheritParams revert
#' @param B number of resamples.
#' @param probs quantiles for the interval.
#' @param seed optional random seed.
#' @return a data frame, one row per landmark in the order of
#'     \code{phases$durations}: the estimated \code{age} and
#'     \code{duration} until the next landmark, and their quantiles over
#'     resamples (\code{age.lo}, \code{age.hi}, \code{duration.lo},
#'     \code{duration.hi}); attribute \code{"ages"}: the landmark ages
#'     of all resamples, landmarks x resamples, for intervals between
#'     any two landmarks.
#' @export
duration_ci <- function(phases, B=100, probs=c(.025, .975), seed) {

    if ( !missing(seed) ) set.seed(seed)
    era <- phases$era
    mk <- phases$marks
    rp <- phases$rotation.phase
    d0 <- phases$durations
    to <- match(d0$to, d0$landmark)

    res.b <- replicate(B, {
        idx <- sample(nrow(rp), replace=TRUE)
        sub <- list(rotation.phase=rp[idx,], era=era)
        ph <- rp[[era$phase]][idx]
        o <- state_phase(mk[, idx], ph)[rownames(mk)[1]]
        lm <- phase_landmarks(sub, mk[, idx], between=phases$between,
                              phase=era$phase)
        d <- phase_durations(sub, lm, origin=o)
        age <- d$age[match(d0$landmark, d$landmark)]
        c(age, ifelse(is.na(to), era$T, age[to]) - age)
    })
    n <- nrow(d0)
    q <- function(x) stats::quantile(x, probs, na.rm=TRUE)
    qa <- apply(res.b[1:n, , drop=FALSE], 1, q)
    qd <- apply(res.b[n + 1:n, , drop=FALSE], 1, q)
    res <- data.frame(landmark=d0$landmark, to=d0$to,
                      age=d0$age, age.lo=qa[1,], age.hi=qa[2,],
                      duration=d0$duration, duration.lo=qd[1,],
                      duration.hi=qd[2,])
    ages <- matrix(res.b[1:n,], nrow=n, dimnames=list(d0$landmark, NULL))
    attr(res, "ages") <- ages
    res
}

#' Time weights of cells along the cycle (ERA).
#'
#' In an exponentially growing population the cells' ages have the
#' density \code{2*log(2)/T * 2^(-age/T)}: young cells are
#' over-represented in a snapshot. Weighting each cell by
#' \code{2^(age/T)}, the inverse of that density, turns averages over
#' cells into averages over the time of one cycle, e.g. the share of
#' its cycle a cell spends in a state.
#' @param age cell ages.
#' @param T doubling time, in the units of \code{age}.
#' @return weights, one per cell, summing to 1.
#' @export
time_weights <- function(age, T=1) {
    w <- 2^(age/T)
    w/sum(w)
}

#' Average over the cycle time of a per-cell value.
#'
#' The mean of \code{y} weighted by \link{time_weights}: the average a
#' single cell sees over one cycle, as opposed to the snapshot mean over
#' cells.
#' @inheritParams time_weights
#' @param y a per-cell value, e.g. a cohort state, or a logical for the
#'     share of time in a state.
#' @return the time average.
#' @export
time_mean <- function(y, age, T=1) {
    ok <- !is.na(y) & !is.na(age)
    sum(time_weights(age[ok], T)*y[ok])
}
