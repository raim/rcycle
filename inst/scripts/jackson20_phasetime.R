### DIVISION-CYCLE PSEUDOPHASE IN TIME, JACKSON ET AL. 2020

## Cell age by ergodic rate analysis (ERA), landmark and group
## durations, per-cell quality, for the 11 media of Jackson et al. 2020.
## Input: a light extract of the owner's pseudophase pipeline, list J
## with per medium counts, total, states, cohorts, cdc.srt and mu
## (growth rate, natural log per hour).
## Usage: Rscript jackson20_phasetime.R [jackson20.rda]
## Output: the list R, saved to jackson20_phasetime.rda.

library(rcycle)
options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly=TRUE)
infile <- "~/work/teaching/pwm/scmw/jackson20.rda"
if ( length(args) ) infile <- args[1]
load(infile)

## division-cycle cohorts in cycle order, and the G1/S point of
## maximal change between the G1 and S peaks
cdc <- c("cdc_M/G1", "cdc_G1", "cdc_S", "cdc_G2", "cdc_M")
between <- list(G1S=c("cdc_G1", "cdc_S"))

## pseudophase with origin at the M/G1 peak (division), S before G2
cycle_phase <- function(st) {
    P <- get_pseudophase(st, revert.phase=TRUE, verb=0)
    pk <- state_phase(st, P$rotation.phase$phi)
    if ( phase_since(pk["cdc_G2"], pk["cdc_M/G1"]) <
         phase_since(pk["cdc_S"], pk["cdc_M/G1"]) )
        P <- revert(P, force=TRUE, verb=0)
    P
}
phase_since <- function(phi, origin) (phi - origin) %% (2*pi)

one_medium <- function(x, B=20, boot=0) {

    st <- x$states[cdc,]
    n <- lengths(x$cohorts[cdc])
    T <- ifelse(x$mu > 0, log(2)/x$mu, NA)

    P <- cycle_phase(st)
    P <- get_segments(P, method='inflection', names='slope', spar=.001,
                      plot=FALSE)
    P <- phase_time(P, T=ifelse(is.na(T), 1, T), origin="cdc_M/G1",
                    states=st)
    lm <- phase_landmarks(P, st, between=between)
    dur <- phase_durations(P, lm)

    ## are the cohort peaks in cycle order?
    order.ok <- identical(dur$landmark[dur$landmark %in% cdc], cdc)

    ## quality against Poisson noise of the same totals
    P <- cell_quality(P, st, x$total, n, B=B, seed=1)

    ## the same, without noise cells
    ok <- !P$rotation.phase$noise
    Pc <- cycle_phase(st[, ok])
    Pc <- get_segments(Pc, method='inflection', names='slope', spar=.001,
                       plot=FALSE)
    Pc <- phase_time(Pc, T=P$era$T, origin="cdc_M/G1", states=st[, ok])
    durc <- phase_durations(Pc, phase_landmarks(Pc, st[, ok],
                                                 between=between))

    ## daughters by their daughter-specific gene counts
    dcnt <- round(x$states["daughter",]*length(x$cohorts$daughter)*x$total)
    sp <- split_counts(dcnt, x$total)
    grp <- ifelse(sp$high, "daughter", "mother")
    gd <- lapply(c(S="cdc_S", G1S="G1S"), function(l)
        group_durations(P, group=grp, to=lm[l]))

    ## per-cell phase uncertainty by resampling genes
    if ( boot > 0 )
        P <- phase_bootstrap(P, x$counts, x$cohorts[cdc], total=x$total,
                             B=boot, seed=2)

    list(T=T, mu=x$mu, n=ncol(st), P=P, landmarks=lm, durations=dur,
         order.ok=order.ok,
         durations.clean=durc, daughter=sp, groups=gd,
         noise=mean(P$rotation.phase$noise))
}

if ( !exists("boot.medium") ) boot.medium <- "jackson20_YPD"
R <- list()
for ( m in names(J) ) {
    cat(paste("medium", m, "\n"))
    R[[sub("jackson20_", "", m)]] <-
        one_medium(J[[m]], boot=ifelse(m %in% boot.medium, 20, 0))
}

## summary table per medium
age_at <- function(d, l) d$age[d$landmark==l]
tab <- do.call(rbind, lapply(names(R), function(m) {
    r <- R[[m]]
    data.frame(medium=m, mu=r$mu, T=r$T, cells=r$n, noise=r$noise,
               order.ok=r$order.ok,
               G1.S=age_at(r$durations, "cdc_S"),
               G1.G1S=age_at(r$durations, "G1S"),
               G1.S.clean=age_at(r$durations.clean, "cdc_S"),
               G1.G1S.clean=age_at(r$durations.clean, "G1S"),
               daughters=mean(r$daughter$high),
               G1d=r$groups$G1S$duration[1], G1m=r$groups$G1S$duration[2],
               G1d.S=r$groups$S$duration[1], G1m.S=r$groups$S$duration[2])
}))
## durations in units of T where T is unknown (no growth)
tab[is.na(tab$T), grep("^G1", colnames(tab))] <- NA
if ( !exists("outfile") ) outfile <- "jackson20_phasetime.rda"
save(R, tab, file=outfile)
print(tab, digits=2)
