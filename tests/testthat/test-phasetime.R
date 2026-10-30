## phase_time and friends (2026-10-30): cell age by ergodic rate analysis
## (ERA, Kafri 2013) from the cell order, landmark and group durations,
## per-cell quality against Poisson noise, closed principal curve.
## Simulated asynchronous exponentially growing populations.

library(testthat)
library(rcycle)

## ages of an asynchronous exponential population, density 2^(-a/T)
sim_ages <- function(n, T=1) -T*log2(1 - runif(n)/2)

## cohort states of cells at ages a: cosine profiles peaking at ages
## pk (fractions of T), Poisson counts from totals tot, n genes each
sim_states <- function(a, T=1, pk=(0:4)/5, amp=.8, tot, n=50, p=.01) {
    w <- 2*pi*a/T
    lam <- sapply(w, function(x) p*(1 + amp*cos(x - 2*pi*pk)))
    cnt <- matrix(rpois(length(lam), lam*rep(tot, each=length(pk))),
                  nrow=length(pk))
    st <- cnt/n/rep(tot, each=length(pk))
    dimnames(st) <- list(paste0("c", seq_along(pk)), paste0("x", seq_along(a)))
    st
}

## circular error in units of T
age_err <- function(x, y, T=1) abs((x - y + T/2) %% T - T/2)

test_that("era_age, era_fraction and era_duration are consistent", {
    F <- seq(0, 1, .1)
    expect_equal(era_fraction(era_age(F, T=3), T=3), F)
    expect_equal(era_age(c(0, 1), T=2), c(0, 2))
    expect_equal(era_duration(F, T=2, share=1), era_age(F, T=2))
    ## half the newborn cells, half of all cells: the whole cycle
    expect_equal(era_duration(.5, T=2, share=.5), 2)
})

test_that("ERA recovers ages from the order of exponentially growing cells", {
    set.seed(1)
    a <- sim_ages(5000, T=2)
    F <- (rank(a) - .5)/length(a)
    expect_lt(max(abs(era_age(F, T=2) - a)), .05)
    ## the linear calibration is biased
    expect_gt(max(abs(2*F - a)), .15)
})

set.seed(2)
T <- 1.5
a <- sim_ages(3000, T=T)
tot <- round(exp(rnorm(length(a), log(3000), .3)))
st <- sim_states(a, T=T, tot=tot)
P <- get_pseudophase(st, revert.phase=TRUE, verb=0)

test_that("phase_time recovers ages from the cohort pseudophase", {

    P <- phase_time(P, T=T, origin="c1", states=st)
    err <- age_err(P$rotation.phase$age, a, T)
    expect_lt(median(err), .03*T)
    expect_equal(P$era$T, T)
    expect_true(all(c("age", "age_F") %in% colnames(P$rotation.phase)))
    expect_true("age" %in% colnames(P$x.phase))

    ## about the same age from theta, whose order is the same
    Pt <- phase_time(P, T=T, origin="c1", phase="theta", states=st)
    expect_lt(median(age_err(Pt$rotation.phase$age,
                             P$rotation.phase$age, T)), .01*T)

    ## landmarks: cohort peaks at known ages
    lm <- phase_landmarks(P, st)
    dur <- phase_durations(P, lm)
    expect_equal(dur$landmark, paste0("c", 1:5))
    expect_lt(max(age_err(dur$age, T*(0:4)/5, T)), .05*T)
    expect_equal(sum(dur$duration), T)
    expect_equal(dur$to, c(paste0("c", 2:5), "division"))
})

test_that("phase_landmarks finds a point of maximal change between peaks", {
    lm <- phase_landmarks(P, st, between=list(c12=c("c1", "c2")))
    expect_equal(names(lm), c(paste0("c", 1:5), "c12"))
    d <- phase_since(lm["c12"], lm["c1"])
    expect_true(d >= 0 & d <= phase_since(lm["c2"], lm["c1"]))
})

## agent-based budding yeast: daughters and mothers both start in G1,
## daughters stay longer; S/G2/M takes B for both; durations vary
## lognormally with coefficient of variation cv
sim_yeast <- function(n0=200, g1=c(d=.6, m=.2), B=1, cv=.1, tend=9) {
    rl <- function(n, m) m*exp(rnorm(n, -cv^2/2, cv))
    birth <- runif(n0, -2, 0)
    type <- sample(c("d", "m"), n0, replace=TRUE)
    out <- list()
    while ( length(birth) ) {
        tg <- rl(length(birth), g1[type])
        div <- birth + tg + rl(length(birth), B)
        out[[length(out)+1]] <- data.frame(birth=birth, type=type,
                                            g1=tg, div=div)
        b <- div[div <= tend]
        birth <- c(b, b)
        type <- rep(c("m", "d"), each=length(b))
    }
    do.call(rbind, out)
}

test_that("group_durations recovers daughter and mother G1 durations", {
    set.seed(7)
    g1 <- c(d=.6, m=.2)
    cv <- .1
    tend <- 9
    all <- sim_yeast(g1=g1, cv=cv, tend=tend)
    alive <- function(t) all$birth <= t & all$div > t
    cells <- all[alive(tend),]
    r <- log(sum(alive(tend))/sum(alive(tend-2)))/2
    Td <- log(2)/r

    ## expected: -log(E[exp(-r*tau)])/r, tau the G1 durations
    x <- exp(rnorm(1e5, -cv^2/2, cv))
    expd <- sapply(g1, function(m) -log(mean(exp(-r*m*x)))/r)

    ## a minimal phases object: the age is the phase, origin 0
    age <- tend - cells$birth
    ph <- list(rotation.phase=data.frame(phi=2*pi*age/(1.01*max(age))),
               era=list(T=Td, origin=0, phase="phi"))
    grp <- ifelse(age < cells$g1, cells$type, NA)
    gd <- group_durations(ph, group=grp, to=2*pi*.999)
    expect_equal(gd$group, c("d", "m", "all"))
    expect_equal(gd$share[1:2], c(.5, .5))
    expect_equal(gd$duration[1:2], unname(expd), tolerance=.05)
})

test_that("split_counts separates two Poisson rates", {
    set.seed(3)
    tot <- round(runif(4000, 1000, 5000))
    hi <- runif(4000) < .3
    x <- rpois(4000, ifelse(hi, 2e-3, 1e-4)*tot)
    sp <- split_counts(x, tot)
    expect_equal(sp$rate, c(1e-4, 2e-3), tolerance=.1)
    expect_gt(mean(sp$high == hi), .95)
})

test_that("project_cells reproduces the PCA loadings", {
    pr <- project_cells(P, st)
    expect_equal(unname(pr$PC1), unname(P$rotation[,1]), tolerance=1e-8)
    expect_equal(unname(pr$PC2), unname(P$rotation[,2]), tolerance=1e-8)
    expect_true(all(pr$r2 >= 0 & pr$r2 <= 1 + 1e-12))
    expect_equal(pr$sd, unname(P$scale), tolerance=1e-8)
})

test_that("cell_quality flags noise cells but not cycling cells", {
    ## add noise cells: same totals, no cycle
    n <- setNames(rep(50, 5), rownames(st))
    nz <- noise_states(st[, 1:500], tot[1:500], n, B=1)
    colnames(nz) <- paste0("noise", 1:500)
    st2 <- cbind(st, nz)
    tot2 <- c(tot, tot[1:500])
    Q <- get_pseudophase(st2, verb=0)
    Q <- cell_quality(Q, st2, tot2, n, B=20, seed=4)
    isn <- grepl("^noise", colnames(st2))
    expect_gt(mean(Q$rotation.phase$noise[isn]), .8)
    expect_lt(mean(Q$rotation.phase$noise[!isn]), .1)
    expect_equal(nrow(Q$noise), 20*ncol(st2))
})

test_that("phase_bootstrap gives smaller phase spread for more counts", {
    set.seed(5)
    ng <- 20
    a <- sim_ages(400)
    w <- 2*pi*a
    genes <- paste0("g", 1:(5*ng))
    pk <- rep((0:4)/5, each=ng)
    tot <- rep(c(2000, 20000), each=200)
    lam <- sapply(seq_along(w), function(j)
        2e-3*(1 + .8*cos(w[j] - 2*pi*pk))*tot[j])
    cnt <- matrix(rpois(length(lam), lam), nrow=length(genes),
                  dimnames=list(genes, paste0("x", 1:400)))
    coh <- split(genes, rep(paste0("c", 1:5), each=ng))
    st <- t(sapply(coh, function(g) colMeans(t(t(cnt[g,])/tot))))
    P <- get_pseudophase(st, verb=0)
    P <- phase_time(P, origin="c1", states=st)
    P <- phase_bootstrap(P, cnt, coh, total=tot, B=30, seed=6)
    rp <- P$rotation.phase
    expect_true(all(c("theta.sd", "F.sd", "age.sd") %in% colnames(rp)))
    expect_lt(median(rp$F.sd[tot > 1e4]), median(rp$F.sd[tot < 1e4]))
})

test_that("principal_phase follows a curved closed path", {
    pc <- principal_phase(P, k=3, K=4)
    lam <- pc$rotation.phase$lambda
    th <- pc$rotation.phase$theta
    ## same direction as theta, close to it
    expect_gt(abs(mean(exp(1i*(lam - th)))), .9)
    expect_true(all(lam >= -pi & lam < pi))
    ## ages from lambda are as good as from phi
    pl <- phase_time(pc, T=T, origin="c1", phase="lambda", states=st)
    err <- age_err(pl$rotation.phase$age, a, T)
    expect_lt(median(err), .03*T)
    ## revert also flips lambda
    rv <- revert(pc, force=TRUE, verb=0)
    expect_equal(rv$rotation.phase$lambda, -lam)
})

test_that("the phase-time plots draw without warnings", {
    n <- setNames(rep(50, 5), rownames(st))
    Q <- phase_time(P, T=T, origin="c1", states=st)
    Q <- cell_quality(Q, st, tot, n, B=20, seed=8)
    dur <- phase_durations(Q, phase_landmarks(Q, st))
    pdf(NULL); on.exit(dev.off())
    expect_no_warning(plotPhasePlane(Q))
    expect_no_warning(plotPhasePlane(Q, raw=TRUE, example=NA))
    expect_no_warning(plotPhaseOrder(Q))
    expect_no_warning(plotEra(Q))
    expect_no_warning(plotAgeStates(Q, st, dur))
    ## mirrored PC2 after revert: cells keep their angle in the plot
    xy <- plotPhasePlane(revert(Q, force=TRUE, verb=0))
    rv <- revert(Q, force=TRUE, verb=0)$rotation.phase$theta
    expect_gt(mean(cos(atan2(xy$y, xy$x) - rv)), .99)
})
