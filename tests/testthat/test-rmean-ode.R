## Test the analytic cycle means against numerical integration of the
## pulse width-modulated ODEs (pwmode_*), with growth dilution mu > 0.
##
## * get_rmean_exact: the periodic mean of the ODEs as coded, where mu dilutes
##   in both phases;
## * get_rmean: for the models with phase-switched degradation, the exact mean
##   of a model WITHOUT dilution in the ON phase (gamma = dr + mu only in OFF);
##   it agrees with the coded ODEs only for mu = 0.

library(testthat)
library(rcycle)

## numerical cycle mean: integrate from R = 0 until periodic, then average over
## the last full cycles (trapezoid rule); the phase switches lie on the grid
ode_mean <- function(func, parms, phi, tau, nstep = 1000, ncyc = NULL) {
    loss <- with(as.list(parms), mu*phi*tau + (dr + mu)*(1 - phi)*tau)
    if ( is.null(ncyc) ) ncyc <- max(20, ceiling(40/loss))
    hocf <- function(t) as.numeric(((t/tau) %% 1) < phi - 1e-9)
    times <- seq(0, ncyc*tau, length.out = ncyc*nstep + 1)
    out <- deSolve::ode(y = c(R = 0), times = times, func = func, parms = parms,
                        hocf = hocf, method = "rk4")
    last <- times >= (ncyc - 5)*tau
    t <- out[last, "time"]; r <- out[last, "R"]
    sum(diff(t) * (r[-1] + r[-length(r)])/2) / (max(t) - min(t))
}

## ODE without dilution in the ON phase: the model get_rmean solves exactly
pwmode_k_dr_k0_noondil <- function(time, state, parameters, hocf) {
    kappa <- hocf(time)
    with(as.list(c(state, parameters)), {
        dR = kappa*k + (1 - kappa)*k0 - (1 - kappa)*(dr + mu)*R
        return(list(c(dR)))
    })
}

PARS <- list(c(k = 10, k0 = 2, dr = 1, mu = 0.1, phi = 0.3, tau = 4),
             c(k = 10, k0 = 5, dr = 3, mu = 0.2, phi = 0.6, tau = 1),
             c(k = 5, k0 = 20, dr = 0.3, mu = 0.05, phi = 0.2, tau = 8),
             c(k = 264, k0 = 26.4, dr = 1.7, mu = 0.1, phi = 0.5, tau = 2))
FUNS <- list(k = pwmode_k, dr = pwmode_dr, k_dr = pwmode_k_dr, k_dr_k0 = pwmode_k_dr_k0)

test_that("get_rmean_exact agrees with the ODEs for mu > 0", {
    skip_if_not_installed("deSolve")
    for ( p in PARS ) for ( mod in names(FUNS) ) {
        parms <- p[c("k", "k0", "dr", "mu")]
        num <- ode_mean(FUNS[[mod]], parms, p[["phi"]], p[["tau"]])
        ana <- get_rmean_exact(k = p[["k"]], k0 = p[["k0"]], dr = p[["dr"]], mu = p[["mu"]],
                               phi = p[["phi"]], tau = p[["tau"]], model = mod)
        expect_equal(ana, num, tolerance = 1e-3,
                     label = paste(mod, paste(names(p), p, sep = "=", collapse = " ")))
    }
})

test_that("get_rmean is exact for the model without dilution in the ON phase", {
    skip_if_not_installed("deSolve")
    for ( p in PARS ) {
        parms <- p[c("k", "k0", "dr", "mu")]
        num <- ode_mean(pwmode_k_dr_k0_noondil, parms, p[["phi"]], p[["tau"]])
        ana <- get_rmean(k = p[["k"]], k0 = p[["k0"]], dr = p[["dr"]], mu = p[["mu"]],
                         phi = p[["phi"]], tau = p[["tau"]], model = "k_dr_k0")
        expect_equal(ana, num, tolerance = 1e-3)
    }
})

test_that("get_rmean_exact and get_rmean agree for mu = 0 and for model k", {
    taus <- seq(0.5, 8, 0.5)
    for ( mod in c("k", "dr", "k_dr", "k_dr_k0") ) {
        x <- get_rmean_exact(k = 264, k0 = 26.4, dr = 1.7, mu = 0, phi = 0.5, tau = taus, model = mod)
        y <- get_rmean(k = 264, k0 = 26.4, dr = 1.7, mu = 0, phi = 0.5, tau = taus, model = mod)
        expect_equal(x, y, tolerance = 1e-8)
    }
    x <- get_rmean_exact(k = 10, dr = 1, mu = 0.2, phi = 0.3, tau = 4, model = "k")
    expect_equal(x, 0.3*10/1.2)
})

test_that("get_rmean_exact is stable for very small mu", {
    x <- get_rmean_exact(k = 10, k0 = 2, dr = 1, mu = c(0, 1e-12, 1e-9, 1e-6), phi = 0.3, tau = 4,
                         model = "k_dr_k0")
    expect_true(all(is.finite(x)))
    expect_equal(x[-1], rep(x[1], 3), tolerance = 1e-5)
})
