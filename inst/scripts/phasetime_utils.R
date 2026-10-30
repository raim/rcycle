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

## pseudophase with origin at the M/G1 peak (division), S before G2
cycle_phase <- function(st, mk) {
    P <- get_pseudophase(st, revert.phase=TRUE, verb=0)
    pk <- state_phase(mk, P$rotation.phase$phi)
    if ( phase_since(pk["cdc_G2"], pk["cdc_M/G1"]) <
         phase_since(pk["cdc_S"], pk["cdc_M/G1"]) )
        P <- revert(P, force=TRUE, verb=0)
    P
}
phase_since <- function(phi, origin) (phi - origin) %% (2*pi)

one_medium <- function(x, B=20, boot=0) {

    st <- x$states[phase.cohorts,]
    mk <- marks(x)
    n <- lengths(x$cohorts[phase.cohorts])
    T <- ifelse(isTRUE(x$mu > 0), log(2)/x$mu, NA)

    P <- cycle_phase(st, mk)
    P <- get_segments(P, method='inflection', names='slope', spar=.001,
                      plot=FALSE)
    P <- phase_time(P, T=ifelse(is.na(T), 1, T), origin="cdc_M/G1",
                    states=mk)
    lm <- phase_landmarks(P, mk, between=between)
    dur <- phase_durations(P, lm)

    ## are the cohort peaks in cycle order?
    order.ok <- identical(dur$landmark[dur$landmark %in% cdc], cdc)

    ## quality against Poisson noise of the same totals
    P <- cell_quality(P, st, x$total, n, B=B, seed=1)

    ## the same, without noise cells
    ok <- !P$rotation.phase$noise
    Pc <- cycle_phase(st[, ok], mk[, ok])
    Pc <- get_segments(Pc, method='inflection', names='slope', spar=.001,
                       plot=FALSE)
    Pc <- phase_time(Pc, T=P$era$T, origin="cdc_M/G1", states=mk[, ok])
    durc <- phase_durations(Pc, phase_landmarks(Pc, mk[, ok],
                                                 between=between))
    width <- sapply(cdc, function(k)
        peak_width(Pc$rotation.phase$age, mk[k, ok]))

    ## daughters by their daughter-specific gene counts
    dcnt <- round(x$states["daughter",]*length(x$cohorts$daughter)*x$total)
    sp <- split_counts(dcnt, x$total)
    grp <- ifelse(sp$high, "daughter", "mother")
    gd <- lapply(c(S="cdc_S", G1S="G1S"), function(l)
        group_durations(P, group=grp, to=lm[l]))

    ## per-cell phase uncertainty by resampling genes
    if ( boot > 0 )
        P <- phase_bootstrap(P, x$counts, x$cohorts[phase.cohorts],
                             total=x$total, B=boot, seed=2)

    list(T=T, mu=x$mu, n=ncol(st), P=P, landmarks=lm, durations=dur,
         order.ok=order.ok,
         durations.clean=durc, width.clean=width, daughter=sp, groups=gd,
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
               ## phases between points of maximal change, and from
               ## peak to peak, without noise cells
               S.clean=dt(r$durations.clean, "G1S", "SG2"),
               G2.clean=dt(r$durations.clean, "SG2", "G2M"),
               SG2pk.clean=dt(r$durations.clean, "cdc_S", "cdc_G2"),
               G2Mpk.clean=dt(r$durations.clean, "cdc_G2", "cdc_M"),
               wS.clean=unname(r$width.clean["cdc_S"]),
               wG2.clean=unname(r$width.clean["cdc_G2"]),
               daughters=mean(r$daughter$high),
               G1d=r$groups$G1S$duration[1], G1m=r$groups$G1S$duration[2],
               G1d.S=r$groups$S$duration[1], G1m.S=r$groups$S$duration[2])
}))
