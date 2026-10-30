### PLOTS: FROM COHORT STATES TO CELL AGE

## Four steps, one plot each: the cells in the PC1/PC2 plane and their
## angle (plotPhasePlane), the angle against the rank phase
## (plotPhaseOrder), the ERA transform of the rank fraction into age
## (plotEra), and the cohort states along age with landmarks
## (plotAgeStates).

## age or phase colors
phase_colors <- function(x, n=100, lim=range(x, na.rm=TRUE))
    grDevices::hcl.colors(n, "viridis")[
                   pmin(n, 1 + floor((n-1)*(x - lim[1])/diff(lim)))]

#' Cells and cohorts in the PC1/PC2 plane, and the phase angle.
#'
#' Cells are the PCA loadings (\code{phases$rotation}), colored by a
#' column of \code{rotation.phase}, e.g. the age from \link{phase_time};
#' the cohorts are drawn as arrows in the direction of their scores
#' (\code{phases$x}). With \code{raw=TRUE}, each cell is drawn at its
#' unscaled amplitude (\link{cell_quality}), so cells without a cycle
#' fall to the center. The angle \code{theta} of one cell is marked.
#' @inheritParams revert
#' @param color column of \code{rotation.phase} for the cell colors.
#' @param raw draw cells at their unscaled amplitude
#'     (\code{amplitude.raw}, requires \link{cell_quality}).
#' @param noise draw cells flagged as noise (\link{cell_quality}) as
#'     grey open circles.
#' @param example index of the cell whose angle is marked; \code{NA}:
#'     none.
#' @param labels cohort labels; default: the row names of
#'     \code{phases$x} without prefix.
#' @param cex point size.
#' @param ... further arguments to \code{plot}.
#' @export
plotPhasePlane <- function(phases, color='age', raw=FALSE, noise=TRUE,
                           example, labels, cex=.3, ...) {

    rp <- phases$rotation.phase
    X <- phases$rotation[,1]
    Y <- phases$rotation[,2]
    ## after revert, theta is mirrored: mirror PC2 likewise
    s <- sign(mean(cos(rp$theta - atan2(Y, X))) -
              mean(cos(rp$theta - atan2(-Y, X))))
    if ( s==0 ) s <- 1
    Y <- s*Y
    if ( raw ) {
        r <- rp$amplitude.raw/rp$amplitude
        X <- X*r
        Y <- Y*r
    }
    col <- phase_colors(rp[[color]])
    nz <- rep(FALSE, length(X))
    if ( noise & 'noise' %in% colnames(rp) ) nz <- rp$noise

    lim <- c(-1, 1)*max(abs(c(X, Y)))
    plot(1, type='n', asp=1, xlim=lim, ylim=lim,
         xlab="PC1 loading", ylab="PC2 loading", ...)
    ## noise cells first, below the others
    if ( any(nz) )
        points(X[nz], Y[nz], col="#00000033", pch=1, cex=cex)
    points(X[!nz], Y[!nz], col=col[!nz], pch=19, cex=cex)
    abline(h=0, v=0, col="grey", lwd=.5)

    ## cohorts: directions of their scores, same mirror as the cells
    cX <- phases$x[,1]
    cY <- s*phases$x[,2]
    ca <- atan2(cY, cX)
    if ( missing(labels) ) labels <- sub(".*_", "", rownames(phases$x))
    L <- .95*lim[2]
    arrows(0, 0, L*cos(ca), L*sin(ca), length=.05, lwd=1.5)
    text(1.08*L*cos(ca), 1.08*L*sin(ca), labels, xpd=TRUE, cex=.8)

    ## the angle of one cell
    if ( missing(example) )
        example <- which.max(rp$amplitude*(abs(rp$theta - 2) < .3))
    if ( !is.na(example) ) {
        th <- atan2(Y[example], X[example])
        segments(0, 0, X[example], Y[example], col=2, lwd=1.5)
        a <- seq(0, th, length.out=30)
        r <- .3*sqrt(X[example]^2 + Y[example]^2)
        lines(r*cos(a), r*sin(a), col=2)
        text(1.4*r*cos(th/2), 1.4*r*sin(th/2), expression(theta), col=2)
        points(X[example], Y[example], col=2, pch=19, cex=4*cex)
    }
    invisible(list(x=X, y=Y))
}

#' Phase angle against rank phase.
#'
#' Each cell's angle \code{theta} against its rank phase \code{phi},
#' the order of the cells along the cycle. A steep stretch is a region
#' of the circle the cells pass quickly (few cells per angle), a flat
#' one a region where cells accumulate. Cohort phases on the top axis.
#' @inheritParams plotPhasePlane
#' @param ... further arguments to \code{plot}.
#' @export
plotPhaseOrder <- function(phases, color='age', cex=.3, ...) {

    rp <- phases$rotation.phase
    o <- order(rp$phi)
    th <- remove_jumps(rp$theta[o], center=FALSE, verb=0)
    plot(rp$phi[o], th, col=phase_colors(rp[[color]])[o], pch=19,
         cex=cex, axes=FALSE,
         xlab=expression(rank~phase~phi), ylab=expression(angle~theta),
         ...)
    circ.axis(1)
    axis(2)
    box()
    axis(3, at=phases$x.phase$phi,
         labels=sub(".*_", "", rownames(phases$x.phase)), las=2,
         cex.axis=.7)
}

#' The ERA transform: fraction of cells to age.
#'
#' The age of a cell with a fraction \code{F} of all cells between
#' division and itself, \code{T*(1-log2(2-F))} (\link{era_age}), against
#' the linear calibration \code{F*T}; the histogram of ages shows the
#' exponentially falling age distribution, density \code{2^(-age/T)}.
#' @inheritParams revert
#' @param T doubling time; default: from \link{phase_time}.
#' @param ... further arguments to \code{plot}.
#' @export
plotEra <- function(phases, T=phases$era$T, ...) {

    era <- phases$era
    F <- phases$rotation.phase[[paste0(era$name, '_F')]]
    age <- phases$rotation.phase[[era$name]]

    h <- graphics::hist(age, breaks=20, plot=FALSE)
    plot(1, type='n', xlim=c(0, 1), ylim=c(0, T),
         xlab="fraction of cells since division, F",
         ylab="age", ...)
    ## age histogram on the right, as horizontal bars
    w <- .3*h$density/max(h$density)
    rect(1 - w, h$breaks[-length(h$breaks)], 1, h$breaks[-1],
         col="#00000020", border=NA)
    Fx <- seq(0, 1, length.out=200)
    lines(Fx, era_age(Fx, T), lwd=2, col=2)
    abline(0, T, col="grey", lty=2)
    legend("topleft", c("ERA", "linear", "density of ages"),
           col=c(2, "grey", "#00000040"), lty=c(1, 2, NA),
           pch=c(NA, NA, 15), pt.cex=1.5, lwd=c(2, 1, NA), bty='n',
           cex=.8)
}

#' Cohort states along age, with landmarks and phases.
#'
#' Moving averages of the cohort states over cells sorted by age (or
#' any phase column), each scaled to its maximum; landmarks
#' (\link{phase_durations}) on the top axis and the phases between
#' them shaded.
#' @inheritParams phase_time
#' @param states cohort states, cohorts x cells.
#' @param durations a table from \link{phase_durations}.
#' @param x column of \code{rotation.phase} for the x-axis.
#' @param window fraction of cells for the moving average.
#' @param col cohort colors.
#' @param shade landmark names that start shaded phases, alternating.
#' @param xlab x-axis label.
#' @param ... further arguments to \code{matplot}.
#' @export
plotAgeStates <- function(phases, states, durations, x='age',
                          window=.05, col=seq_len(nrow(states)),
                          shade, xlab=x, ...) {

    a <- phases$rotation.phase[[x]]
    o <- order(a)
    n <- max(1, round(window*length(a)))
    ma <- t(apply(states[, o, drop=FALSE], 1, function(y)
        stats::filter(y, rep(1, n)/n, circular=TRUE)))
    ma <- ma/apply(ma, 1, max)

    matplot(a[o], t(ma), type='n', xlab=xlab, ylab="state / max", ...)
    if ( !missing(durations) ) {
        d <- durations
        if ( missing(shade) ) shade <- d$landmark
        for ( i in which(d$landmark %in% shade)[c(TRUE, FALSE)] )
            rect(d$age[i], 0, d$age[i] + d$duration[i], 2,
                 col="#00000012", border=NA)
        abline(v=d$age, col="grey", lwd=.5)
        axis(3, at=d$age, labels=sub("cdc_", "", d$landmark), las=2,
             cex.axis=.7)
    }
    matlines(a[o], t(ma), lty=1, col=col)
    legend("bottomleft", sub(".*_", "", rownames(states)), col=col,
           lty=1, cex=.7, bty='n')
}
