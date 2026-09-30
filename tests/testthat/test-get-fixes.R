## Tests for four fixes in the get_* functions (2026-09-30):
## 1. get_basal/get_rates: an Rmax (or Rmin) that is NA counts as not given;
##    get_rates returned k0 = NA whenever only Rmin was given
## 2. get_ramp: gamma from dr + mu also when k is not given
## 3. get_pmean: RNA rates are passed to get_rmean via ..., not taken from
##    the calling environment
## 4. get_tau/get_times: the root search for the models with phase-switched
##    degradation starts above the pole tau = A/k (phi = 1)
## Reference values: the analytic solutions at mu = 0, where they are exact
## for the pwmode_* ODEs (see test-rmean-ode.R).

library(testthat)
library(rcycle)

## a periodic K_DR_K0 steady state at mu = 0: mean, min, max, amplitude
k <- 10; k0 <- 2; dr <- 1; phi <- 0.3; tau <- 4
Rmean <- get_rmean(k = k, k0 = k0, dr = dr, mu = 0, phi = phi, tau = tau, model = "k_dr_k0")
Rmin <- k*phi*tau/expm1(dr*tau*(1 - phi)) + k0/dr
Rmax <- Rmin + k*phi*tau

test_that("get_basal returns k0 from Rmin, Rmax or both, ignoring NA", {
    expect_equal(get_basal(k = k, gamma = dr, tau = tau, phi = phi, Rmin = Rmin), k0)
    expect_equal(get_basal(k = k, gamma = dr, tau = tau, phi = phi, Rmax = Rmax), k0)
    expect_equal(get_basal(k = k, gamma = dr, tau = tau, phi = phi, Rmin = Rmin, Rmax = Rmax,
                           verb = 0), k0)
    expect_equal(get_basal(k = k, gamma = dr, tau = tau, phi = phi, Rmin = Rmin, Rmax = NA), k0)
    expect_equal(get_basal(k = k, gamma = dr, tau = tau, phi = phi, Rmin = NA, Rmax = Rmax), k0)
    ## element-wise over vectors with some NA
    expect_equal(get_basal(k = k, gamma = dr, tau = tau, phi = phi,
                           Rmin = c(Rmin, NA, Rmin), Rmax = c(NA, Rmax, Rmax), verb = 0),
                 rep(k0, 3))
})

test_that("get_rates recovers k0 of K_DR_K0 from Rmin alone", {
    r <- get_rates(model = "k_dr_k0", a = (Rmax - Rmin)/Rmean, R = Rmean, Rmin = Rmin,
                   phi = phi, tau = tau, mu = 0)
    expect_equal(r$k, k, tolerance = 1e-6)
    expect_equal(r$dr, dr, tolerance = 1e-6)
    expect_equal(r$k0, k0, tolerance = 1e-6)
})

test_that("get_ramp without k computes gamma from dr and mu", {
    for ( mod in c("k_dr", "dr") ) {
        x <- get_ramp(dr = 1, mu = 0.1, phi = phi, tau = tau, relative = TRUE, model = mod)
        y <- get_ramp(gamma = 1.1, phi = phi, tau = tau, relative = TRUE, model = mod)
        expect_equal(x, y)
    }
})

test_that("get_pmean with R missing uses the rates passed, not the environment", {
    mu <- 0.1; rho <- 5; l <- 2; dp <- 0.05
    R <- get_rmean(k = k, k0 = k0, dr = dr, mu = mu, phi = phi, tau = tau, model = "k_dr_k0")
    p.ref <- get_pmean(R = R, rho = rho, l = l, dp = dp, mu = mu)
    ## decoys in the calling and global environment must not be used
    local({
        dr <- 99; phi <- 0.99
        p <- get_pmean(rho = rho, l = l, dp = dp, mu = mu, r.model = "k_dr_k0",
                       k = 10, k0 = 2, dr = 1, phi = 0.3, tau = 4)
        expect_equal(p, p.ref)
    })
})

test_that("get_times finds tau above the pole tau = A/k", {
    ## short period, long duty cycle: tau = 1, phi = 0.6, gamma = 3
    for ( mod in c("k_dr", "dr", "k_dr_k0") ) {
        kk <- 10; g <- 3; ph <- 0.6; ta <- 1; kk0 <- if ( mod == "dr" ) kk else if ( mod == "k_dr" ) 0 else 5
        R <- get_rmean(k = kk, k0 = kk0, dr = g, mu = 0, phi = ph, tau = ta, model = mod)
        Rmn <- kk*ph*ta/expm1(g*ta*(1 - ph)) + kk0/g
        a <- kk*ph*ta/R
        tm <- get_times(model = mod, a = a, R = R, Rmin = Rmn, k = kk, gamma = g)
        expect_equal(tm$tau, ta, tolerance = 1e-5, label = mod)
        expect_equal(tm$phi, ph, tolerance = 1e-5, label = mod)
    }
})
