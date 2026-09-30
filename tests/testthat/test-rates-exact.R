## Exact periodic steady state, amplitude and rate inversion with dilution in
## both phases (get_rcycle_exact, get_ramp_exact, get_rates_exact), tested
## against numerical integration of the pwmode_* ODEs, and against the
## closed forms (get_ramp, get_rates) at mu = 0, where those are exact.

library(testthat)
library(rcycle)

## periodic steady state by numerical integration: R at the start (R0) and
## end (R1) of the ON phase, min, max and mean over the last cycle
ode_cycle <- function(model, k, k0, dr, mu, phi, tau, nstep = 5000) {
    func <- list(k = pwmode_k, dr = pwmode_dr, k_dr = pwmode_k_dr,
                 k_dr_k0 = pwmode_k_dr_k0)[[model]]
    loss <- mu*phi*tau + (dr + mu)*(1 - phi)*tau
    ncyc <- max(20, ceiling(40/loss))
    hocf <- function(t) as.numeric(((t/tau) %% 1) < phi - 1e-9)
    times <- seq(0, ncyc*tau, length.out = ncyc*nstep + 1)
    out <- deSolve::ode(y = c(R = 0), times = times, func = func,
                        parms = c(k = k, k0 = k0, dr = dr, mu = mu),
                        hocf = hocf, method = "rk4")
    last <- which(times >= (ncyc - 1)*tau - tau/nstep/2)
    t <- out[last, "time"]; r <- out[last, "R"]
    ## R0 at the end of the last cycle (= start of the next ON phase),
    ## R1 at the grid point closest to the end of the last ON phase
    c(R0 = r[length(r)], R1 = r[which.min(abs(t - ((ncyc - 1)*tau + phi*tau)))],
      Rmin = min(r), Rmax = max(r),
      mean = sum(diff(t)*(r[-1] + r[-length(r)])/2)/tau)
}

## k, k0, dr, mu, phi, tau; the last two in the high-basal regime
## (normal: k/mu > k0/gamma; reversed: k0/gamma > k/mu, R falls during ON)
PARS <- list(c(k = 10, k0 = 2, dr = 1, mu = 0.1, phi = 0.3, tau = 4),
             c(k = 10, k0 = 5, dr = 3, mu = 0.2, phi = 0.6, tau = 1),
             c(k = 5, k0 = 20, dr = 0.3, mu = 0.05, phi = 0.2, tau = 8),
             c(k = 1, k0 = 20, dr = 0.3, mu = 0.2, phi = 0.3, tau = 6))

test_that("get_rcycle_exact agrees with the ODEs, incl. the reversed regime", {
    skip_if_not_installed("deSolve")
    for ( p in PARS ) for ( mod in c("k", "dr", "k_dr", "k_dr_k0") ) {
        num <- ode_cycle(mod, p[["k"]], p[["k0"]], p[["dr"]], p[["mu"]], p[["phi"]], p[["tau"]])
        ana <- get_rcycle_exact(k = p[["k"]], k0 = p[["k0"]], dr = p[["dr"]], mu = p[["mu"]],
                                phi = p[["phi"]], tau = p[["tau"]], model = mod)
        lab <- paste(mod, paste(names(p), p, sep = "=", collapse = " "))
        ## the values exactly at the (discontinuous) switches carry an
        ## O(step) error in the numerical solution; min, max and mean do not
        for ( v in c("Rmin", "Rmax", "mean") )
            expect_equal(ana[[v]], num[[v]], tolerance = 1e-3, label = paste(lab, v))
        ## direction: R rises during ON (normal) or falls (reversed), in both
        expect_equal(sign(ana$R1 - ana$R0), sign(num[["R1"]] - num[["R0"]]), label = lab)
    }
    ## the last parameter set is reversed for k_dr_k0: R falls during ON
    p <- PARS[[4]]
    cy <- get_rcycle_exact(k = p[["k"]], k0 = p[["k0"]], dr = p[["dr"]], mu = p[["mu"]],
                           phi = p[["phi"]], tau = p[["tau"]], model = "k_dr_k0")
    expect_gt(cy$R0, cy$R1)
})

test_that("get_ramp_exact agrees with the ODEs and with get_ramp for mu = 0", {
    skip_if_not_installed("deSolve")
    for ( p in PARS[1:3] ) for ( mod in c("k", "dr", "k_dr", "k_dr_k0") ) {
        num <- ode_cycle(mod, p[["k"]], p[["k0"]], p[["dr"]], p[["mu"]], p[["phi"]], p[["tau"]])
        A <- get_ramp_exact(k = p[["k"]], k0 = p[["k0"]], dr = p[["dr"]], mu = p[["mu"]],
                            phi = p[["phi"]], tau = p[["tau"]], relative = FALSE, model = mod)
        a <- get_ramp_exact(k = p[["k"]], k0 = p[["k0"]], dr = p[["dr"]], mu = p[["mu"]],
                            phi = p[["phi"]], tau = p[["tau"]], model = mod)
        expect_equal(A, num[["Rmax"]] - num[["Rmin"]], tolerance = 1e-3)
        expect_equal(a, (num[["Rmax"]] - num[["Rmin"]])/num[["mean"]], tolerance = 1e-3)
    }
    for ( mod in c("dr", "k_dr", "k_dr_k0") ) for ( rel in c(TRUE, FALSE) ) {
        x <- get_ramp_exact(k = 10, k0 = 2, dr = 1, mu = 0, phi = 0.3, tau = c(1, 4, 8),
                            relative = rel, model = mod)
        y <- get_ramp(k = 10, k0 = 2, dr = 1, mu = 0, phi = 0.3, tau = c(1, 4, 8),
                      relative = rel, model = mod)
        expect_equal(x, y, tolerance = 1e-8)
    }
    ## model k: get_ramp is exact for any mu
    expect_equal(get_ramp_exact(k = 10, dr = 1, mu = 0.2, phi = 0.3, tau = 4, model = "k"),
                 get_ramp(k = 10, dr = 1, mu = 0.2, phi = 0.3, tau = 4, model = "k"),
                 tolerance = 1e-8)
})

## The inversion is tested on exact forward values (the forward solution is
## tested against the ODEs above): for k_dr_k0 it is ill-conditioned, a 0.1%
## error in Rmin or a can move k0 by several percent (short periods, high dr)
test_that("get_rates_exact inverts get_rcycle_exact, incl. the reversed regime", {
    for ( i in seq_along(PARS) ) for ( mod in c("k", "dr", "k_dr", "k_dr_k0") ) {
        p <- PARS[[i]]
        cy <- get_rcycle_exact(k = p[["k"]], k0 = p[["k0"]], dr = p[["dr"]], mu = p[["mu"]],
                               phi = p[["phi"]], tau = p[["tau"]], model = mod)
        reg <- if ( mod == "k_dr_k0" && i == 4 ) "reversed" else "normal"
        r <- get_rates_exact(model = mod, a = (cy$Rmax - cy$Rmin)/cy$mean, R = cy$mean,
                             Rmin = cy$Rmin, phi = p[["phi"]], tau = p[["tau"]],
                             mu = p[["mu"]], regime = reg, upper = 1e3)
        lab <- paste(mod, reg, paste(names(p), p, sep = "=", collapse = " "))
        expect_equal(r$k, p[["k"]], tolerance = 1e-5, label = paste(lab, "k"))
        expect_equal(r$dr, p[["dr"]], tolerance = 1e-5, label = paste(lab, "dr"))
        if ( mod == "k_dr_k0" )
            expect_equal(r$k0, p[["k0"]], tolerance = 1e-5, label = paste(lab, "k0"))
    }
})

test_that("get_rates_exact agrees with get_rates for mu = 0, and takes Rmax or A", {
    p <- PARS[[1]]
    for ( mod in c("dr", "k_dr", "k_dr_k0") ) {
        cy <- get_rcycle_exact(k = p[["k"]], k0 = p[["k0"]], dr = p[["dr"]], mu = 0,
                               phi = p[["phi"]], tau = p[["tau"]], model = mod)
        a <- (cy$Rmax - cy$Rmin)/cy$mean
        x <- get_rates_exact(model = mod, a = a, R = cy$mean, Rmin = cy$Rmin,
                             phi = p[["phi"]], tau = p[["tau"]], mu = 0)
        y <- get_rates(model = mod, a = a, R = cy$mean, Rmin = cy$Rmin, Rmax = cy$Rmax,
                       phi = p[["phi"]], tau = p[["tau"]], mu = 0)
        expect_equal(x$k, y$k, tolerance = 1e-6)
        expect_equal(x$dr, y$dr, tolerance = 1e-6)
        if ( mod == "k_dr_k0" ) {
            expect_equal(x$k0, y$k0, tolerance = 1e-6)
            z <- get_rates_exact(model = mod, A = cy$Rmax - cy$Rmin, R = cy$mean,
                                 Rmax = cy$Rmax, phi = p[["phi"]], tau = p[["tau"]], mu = 0)
            expect_equal(z$k0, x$k0, tolerance = 1e-6)
        }
    }
})

test_that("get_rates_exact is vectorised and returns NA without a solution", {
    cy <- get_rcycle_exact(k = 10, dr = c(0.5, 1, 2), mu = 0.1, phi = 0.3, tau = 4, model = "k_dr")
    r <- get_rates_exact(model = "k_dr", a = (cy$Rmax - cy$Rmin)/cy$mean, R = cy$mean,
                         phi = 0.3, tau = 4, mu = 0.1)
    expect_equal(r$dr, c(0.5, 1, 2), tolerance = 1e-6)
    expect_equal(r$k, rep(10, 3), tolerance = 1e-6)
    ## a relative amplitude no dr can produce
    r <- get_rates_exact(model = "k_dr", a = 100, R = 1, phi = 0.3, tau = 4, mu = 0.1)
    expect_true(is.na(r$dr))
})
