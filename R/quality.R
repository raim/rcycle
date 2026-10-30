### PER-CELL PHASE QUALITY

## With scale=TRUE every cell is scaled to unit variance across
## cohorts, so all cells have the same amplitude and noise cells still
## get an angle. Here we recover the unscaled amplitude, the share of a
## cell's variance in the PC1/PC2 plane, and compare both with cells
## that have the same total counts but only Poisson noise.

#' Project cohort state columns into the PCA of a phases object.
#'
#' Applies the preprocessing of \link{get_pseudophase} (log, the row
#' means of the original states, column centering and scaling) to new
#' cells and computes their loadings from the cohort scores. For the
#' cells the PCA was built from, this returns their loadings.
#' @inheritParams revert
#' @param states a cohort state matrix with the same rows as the one
#'     used in \link{get_pseudophase}, cells in columns.
#' @param k number of PCs to return.
#' @return a data frame with one row per cell: the loadings
#'     \code{PC1..PCk}; the angle \code{theta} in PC1/PC2 of the original
#'     PCA, i.e. before \link{revert}; the scaled \code{amplitude};
#'     the cell's standard deviation across cohorts after row centering
#'     (\code{sd}); the share of its variance in the PC1/PC2 plane
#'     (\code{r2}); and the unscaled amplitude in that plane,
#'     \code{amplitude.raw = sd*sqrt(r2)}, in units of the states.
#' @export
project_cells <- function(phases, states, k=3) {

    if ( !'row.means' %in% names(phases) )
        stop("phases object without preprocessing, ",
             "re-run get_pseudophase")
    states <- as.matrix(states)[rownames(phases$x), , drop=FALSE]

    z <- states
    if ( isTRUE(phases$log) ) z <- log10(z+1)
    if ( !is.null(phases$row.means) ) z <- z - phases$row.means

    ## column centering and scaling as in prcomp
    if ( !isFALSE(phases$center) )
        z <- t(t(z) - colMeans(z))
    nc <- nrow(z)
    sdv <- sqrt(colSums(z^2)/(nc-1))
    if ( !isFALSE(phases$scale) )
        z <- t(t(z)/sdv)

    ## loadings: z = x %*% t(rotation), the columns of x are orthogonal
    ## with squared norm (n-1)*sdev^2
    sdev <- phases$sdev
    npc <- sum(sdev > max(sdev)*1e-8)
    x <- phases$x[, 1:npc, drop=FALSE]
    r <- t(t(crossprod(z, x))/((nc-1)*sdev[1:npc]^2))

    ## share of variance in PC1/PC2
    vr <- t(t(r^2)*sdev[1:npc]^2)
    r2 <- rowSums(vr[, 1:2, drop=FALSE])/rowSums(vr)

    res <- data.frame(r[, 1:min(k, npc), drop=FALSE],
                      theta=atan2(r[,2], r[,1]),
                      amplitude=sqrt(r[,1]^2 + r[,2]^2),
                      sd=sdv, r2=r2,
                      amplitude.raw=sdv*sqrt(r2))
    colnames(res)[1:min(k, npc)] <- paste0("PC", 1:min(k, npc))
    rownames(res) <- colnames(states)
    res
}

#' Cohort states of noise cells with given total counts.
#'
#' Draws cohort read counts from a Poisson distribution with the
#' population mean fraction of each cohort times the total counts of
#' each cell, i.e. cells without any cycle, and converts them to states
#' as in \link{get_states} (mean fraction per gene of the cohort).
#' Exact if each cohort state is the mean over its genes of counts
#' divided by the cell total, and cohorts share no genes.
#' @param states the cohort state matrix, cohorts x cells.
#' @param total the total counts of each cell (column of
#'     \code{states}).
#' @param n the number of genes per cohort, named by cohort.
#' @param B number of noise cells per cell.
#' @return a state matrix with \code{B} noise cells per cell, named
#'     \code{<cell>.<b>}; attribute \code{"cell"} holds the index of
#'     the original cell.
#' @export
noise_states <- function(states, total, n, B=10) {

    states <- as.matrix(states)
    n <- n[rownames(states)]

    ## population fraction of each cohort's reads
    cnt <- states * n * rep(total, each=nrow(states))
    p <- rowSums(cnt)/sum(total)

    idx <- rep(seq_len(ncol(states)), B)
    tot <- total[idx]
    sim <- matrix(stats::rpois(length(p)*length(idx), p %o% tot),
                  nrow=length(p))
    sim <- sim/n/rep(tot, each=length(p))
    dimnames(sim) <- list(rownames(states),
                          paste0(colnames(states)[idx], ".",
                                 rep(seq_len(B), each=ncol(states))))
    attr(sim, "cell") <- idx
    sim
}

#' Per-cell phase quality against a Poisson noise reference.
#'
#' Adds the unscaled amplitude and the share of variance in the PC1/PC2
#' plane (\link{project_cells}) to \code{phases$rotation.phase}, and
#' compares both with \code{B} noise cells of the same total counts
#' (\link{noise_states}) projected into the same PCA. A cell whose
#' amplitude is not above its noise cells has a phase angle from noise.
#' @inheritParams noise_states
#' @inheritParams revert
#' @param alpha cells with \code{noise.p > alpha} are flagged
#'     \code{noise}.
#' @param seed optional random seed.
#' @return the phases object with columns \code{sd}, \code{r2},
#'     \code{amplitude.raw}, the noise p-values of the amplitude
#'     (\code{noise.p}) and of \code{r2} (\code{noise.p.r2}), and the flag
#'     \code{noise}, in \code{rotation.phase}; the noise cells'
#'     projections in \code{phases$noise}.
#' @export
cell_quality <- function(phases, states, total, n, B=20, alpha=.05,
                         seed) {

    if ( !missing(seed) ) set.seed(seed)

    obs <- project_cells(phases, states)
    sim <- noise_states(states[rownames(phases$x),], total, n, B=B)
    nul <- project_cells(phases, sim)
    cell <- attr(sim, "cell")

    ## empirical p-values per cell
    pval <- function(o, s)
        (1 + tabulate(cell[s >= o[cell]], nbins=length(o)))/(B+1)
    pa <- pval(obs$amplitude.raw, nul$amplitude.raw)
    pr <- pval(obs$r2, nul$r2)

    rp <- phases$rotation.phase
    for ( id in c("sd", "r2", "amplitude.raw") ) rp[[id]] <- obs[[id]]
    rp$noise.p <- pa
    rp$noise.p.r2 <- pr
    rp$noise <- pa > alpha
    phases$rotation.phase <- rp

    nul$cell <- cell
    phases$noise <- nul
    phases$processing <- c(phases$processing, "quality")
    phases
}

## cohort x gene matrix of mean weights, genes may repeat
cohort_weights <- function(idx, genes) {
    m <- matrix(0, nrow=length(idx), ncol=length(genes),
                dimnames=list(names(idx), genes))
    for ( j in seq_along(idx) ) {
        tb <- table(match(idx[[j]], genes))
        m[j, as.integer(names(tb))] <- tb/length(idx[[j]])
    }
    m
}

#' Per-cell phase uncertainty by resampling genes within cohorts.
#'
#' Resamples the genes of each cohort with replacement, recomputes the
#' cohort states from the counts and projects them into the PCA
#' (\link{project_cells}). The spread of a cell's angle over resamples
#' is converted to the spread of its fraction of cells \code{F}, by its
#' rank among the original cells, and, after \link{phase_time}, to the
#' spread of its age.
#' @inheritParams revert
#' @param counts gene counts, genes x cells (also a sparse
#'     \code{Matrix}), with gene names as row names.
#' @param cohorts list of gene names (or row indices) per cohort, with
#'     the cohort names of the PCA.
#' @param total total counts per cell; default: column sums of
#'     \code{counts}.
#' @param B number of resamples.
#' @param seed optional random seed.
#' @return the phases object with columns \code{theta.sd} (circular
#'     standard deviation, radian), \code{F.sd} and, after
#'     \link{phase_time}, \code{<age>.sd} in \code{rotation.phase}.
#' @export
phase_bootstrap <- function(phases, counts, cohorts, total, B=50, seed) {

    if ( !missing(seed) ) set.seed(seed)
    if ( missing(total) ) total <- Matrix::colSums(counts)

    cohorts <- cohorts[rownames(phases$x)]
    genes <- rownames(counts)
    if ( is.numeric(unlist(cohorts)) )
        cohorts <- lapply(cohorts, function(x) genes[x])
    use <- unique(unlist(cohorts))
    counts <- counts[use, , drop=FALSE]

    get_st <- function(idx) {
        st <- as.matrix(cohort_weights(idx, use) %*% counts)
        t(t(st)/total)
    }

    ## original angles in the frame of the PCA, sorted for ranks
    th0 <- project_cells(phases, get_st(cohorts))$theta
    sth <- sort(th0)
    F0 <- findInterval(th0, sth)/length(th0)

    cs <- sn <- dF2 <- 0
    for ( b in seq_len(B) ) {
        idx <- lapply(cohorts, function(x) sample(x, replace=TRUE))
        th <- project_cells(phases, get_st(idx))$theta
        cs <- cs + cos(th - th0)
        sn <- sn + sin(th - th0)
        dF <- findInterval(th, sth)/length(th0) - F0
        dF <- (dF + .5) %% 1 - .5
        dF2 <- dF2 + dF^2
    }
    R <- sqrt(cs^2 + sn^2)/B

    rp <- phases$rotation.phase
    rp$theta.sd <- sqrt(-2*log(pmin(R, 1)))
    rp$F.sd <- sqrt(dF2/B)
    era <- phases$era
    if ( !is.null(era) ) {
        ## d age / dF = T/((2-F)*log(2))
        F <- rp[[paste0(era$name, '_F')]]
        rp[[paste0(era$name, '.sd')]] <- rp$F.sd*era$T/((2-F)*log(2))
    }
    phases$rotation.phase <- rp
    phases$processing <- c(phases$processing, "bootstrap")
    phases
}
