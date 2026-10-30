### DIVISION-CYCLE PSEUDOPHASE IN TIME: SHARED FUNCTIONS

## Used by jackson20_phasetime.R and nadal25_phasetime.R. A data set x
## is a list with counts (genes x cells), total, states (cohorts x
## cells, incl. the phase, landmark and daughter cohorts), cohorts (gene
## lists) and mu (growth rate, natural log per hour; NA: unknown, ages
## are then in units of the cycle).

## cohorts for the phase; landmark cohorts in cycle order, named as
## cdc_*, and the points of maximal change between neighbouring peaks
## (G1/S, S/G2, G2/M)
if ( !exists("phase.cohorts") ) phase.cohorts <- paste0("cdc_bins_", 1:5)
if ( !exists("mark.cohorts") )
    mark.cohorts <- paste0("sp_", c("M/G1", "G1", "S", "G2", "M"))
cdc <- c("cdc_M/G1", "cdc_G1", "cdc_S", "cdc_G2", "cdc_M")
between <- list(G1S=c("cdc_G1", "cdc_S"), SG2=c("cdc_S", "cdc_G2"),
                G2M=c("cdc_G2", "cdc_M"))

## width (full width at half maximum, units of age) of a cohort's
## moving average over w of the cells along age
peak_width <- function(age, y, w=.05) {
    o <- order(age)
    a <- age[o]
    n <- round(w*length(y))
    ma <- stats::filter(y[o], rep(1, n)/n, circular=TRUE)
    ma <- (ma - min(ma))/(max(ma) - min(ma))
    ## age step to each next cell, the last wraps to the first
    da <- diff(c(a, a[1] + max(a)*length(a)/(length(a) - .5)))
    sum(da[ma >= .5])
}

## landmark cohort states, named as cdc
marks <- function(x) {
    mk <- x$states[mark.cohorts,]
    rownames(mk) <- cdc
    mk
}

phase_since <- function(phi, origin) (phi - origin) %% (2*pi)

## interval between two landmarks and its quantiles over the resamples
## of duration_ci
ci_between <- function(ci, a, b, probs=c(.025, .975)) {
    if ( is.null(ci) ) return(rep(NA, 3))
    ag <- attr(ci, "ages")
    est <- ci$age[ci$landmark==b] - ci$age[ci$landmark==a]
    unname(c(est, stats::quantile(ag[b,] - ag[a,], probs, na.rm=TRUE)))
}

ci_at <- function(ci, l, col)
    if ( is.null(ci) ) NA else ci[[col]][ci$landmark==l]

## one data set: phase, ERA age, landmarks and durations with all cells
## and without noise cells (division_phase), CIs of the latter by
## resampling cells, peak widths, daughter calls, and optionally the
## per-cell phase uncertainty by resampling genes
one_medium <- function(x, B=20, boot=0, ci=100) {

    T <- ifelse(isTRUE(x$mu > 0), log(2)/x$mu, NA)
    Tu <- ifelse(is.na(T), 1, T)
    coh <- x$cohorts[phase.cohorts]
    mks <- setNames(x$cohorts[mark.cohorts], cdc)

    P <- division_phase(x$counts, x$total, coh, mks, T=Tu,
                        between=between, B=B, seed=1)
    Pc <- division_phase(x$counts, x$total, coh, mks, T=Tu,
                         between=between, B=B, noise.rm=TRUE, seed=1)
    cic <- NULL
    if ( ci > 0 ) cic <- duration_ci(Pc, B=ci, seed=3)
    width <- sapply(cdc, function(k)
        peak_width(Pc$rotation.phase$age, Pc$marks[k,]))

    ## daughters by their daughter-specific gene counts
    dg <- intersect(x$cohorts$daughter, rownames(x$counts))
    dcnt <- Matrix::colSums(x$counts[dg, , drop=FALSE])
    sp <- split_counts(dcnt, x$total)
    grp <- ifelse(sp$high, "daughter", "mother")
    gd <- lapply(c(S="cdc_S", G1S="G1S"), function(l)
        group_durations(P, group=grp, to=P$landmarks[l]))

    ## per-cell phase uncertainty by resampling genes
    if ( boot > 0 )
        P <- phase_bootstrap(P, x$counts, coh, total=x$total, B=boot,
                             seed=2)

    list(T=T, mu=x$mu, n=ncol(x$counts), P=P, Pc=Pc,
         landmarks=P$landmarks, durations=P$durations,
         order.ok=P$order.ok, durations.clean=Pc$durations,
         ci.clean=cic, width.clean=width, daughter=sp, groups=gd,
         noise=mean(P$rotation.phase$noise))
}

## summary table: one row per data set; durations in hours, or in units
## of the cycle where the growth rate is unknown
age_at <- function(d, l) d$age[d$landmark==l]
dt <- function(d, a, b) age_at(d, b) - age_at(d, a)
phase_table <- function(R) do.call(rbind, lapply(names(R), function(m) {
    r <- R[[m]]
    data.frame(medium=m, mu=r$mu, T=r$T, cells=r$n, noise=r$noise,
               order.ok=r$order.ok,
               G1.S=age_at(r$durations, "cdc_S"),
               G1.G1S=age_at(r$durations, "G1S"),
               G1.S.clean=age_at(r$durations.clean, "cdc_S"),
               G1.G1S.clean=age_at(r$durations.clean, "G1S"),
               G1.lo=ci_at(r$ci.clean, "G1S", "age.lo"),
               G1.hi=ci_at(r$ci.clean, "G1S", "age.hi"),
               ## phases between points of maximal change, and from
               ## peak to peak, without noise cells
               S.clean=dt(r$durations.clean, "G1S", "SG2"),
               S.lo=ci_between(r$ci.clean, "G1S", "SG2")[2],
               S.hi=ci_between(r$ci.clean, "G1S", "SG2")[3],
               G2.clean=dt(r$durations.clean, "SG2", "G2M"),
               SG2pk.clean=dt(r$durations.clean, "cdc_S", "cdc_G2"),
               G2Mpk.clean=dt(r$durations.clean, "cdc_G2", "cdc_M"),
               wS.clean=unname(r$width.clean["cdc_S"]),
               wG2.clean=unname(r$width.clean["cdc_G2"]),
               daughters=mean(r$daughter$high),
               G1d=r$groups$G1S$duration[1], G1m=r$groups$G1S$duration[2],
               G1d.S=r$groups$S$duration[1], G1m.S=r$groups$S$duration[2])
}))
