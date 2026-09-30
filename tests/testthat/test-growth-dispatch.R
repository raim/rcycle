## No growth (mu missing, NA or 0): the exact functions with dilution in both
## phases (get_rmean, get_ramp, get_rates) equal the closed forms
## (get_*_nogrowth), and get_rates uses them; get_rates_nogrowth with
## mu = NA returns the total loss rate as dr, like mu = 0 (it returned NA,
## since get_degradation subtracted mu = NA).

library(testthat)
library(rcycle)

k <- 10; k0 <- 2; dr <- 1; phi <- 0.3; tau <- 4
MODS <- c("dr", "k_dr", "k_dr_k0")
cyc <- function(mod, mu = 0)
    get_rcycle(k = k, k0 = k0, dr = dr, mu = mu, phi = phi, tau = tau, model = mod)

test_that("get_rates_nogrowth with mu = NA equals mu = 0", {
    for ( mod in c("k", MODS) ) {
        cy <- cyc(mod)
        a <- (cy$Rmax - cy$Rmin)/cy$mean
        x <- get_rates_nogrowth(model = mod, a = a, R = cy$mean, Rmin = cy$Rmin,
                                phi = phi, tau = tau, mu = NA)
        y <- get_rates_nogrowth(model = mod, a = a, R = cy$mean, Rmin = cy$Rmin,
                                phi = phi, tau = tau, mu = 0)
        expect_false(any(is.na(x$dr)), label = mod)
        expect_equal(x, y, tolerance = 1e-10, label = mod)
        expect_equal(x$dr, dr, tolerance = 1e-6, label = mod)
    }
})

test_that("get_rates without growth rate uses the closed forms", {
    for ( mod in c("k", MODS) ) {
        cy <- cyc(mod)
        a <- (cy$Rmax - cy$Rmin)/cy$mean
        ref <- get_rates_nogrowth(model = mod, a = a, R = cy$mean, Rmin = cy$Rmin,
                                  phi = phi, tau = tau, mu = 0)
        expect_equal(get_rates(model = mod, a = a, R = cy$mean, Rmin = cy$Rmin,
                               phi = phi, tau = tau), ref, label = paste(mod, "missing"))
        expect_equal(get_rates(model = mod, a = a, R = cy$mean, Rmin = cy$Rmin,
                               phi = phi, tau = tau, mu = NA), ref, label = paste(mod, "NA"))
        expect_equal(get_rates(model = mod, a = a, R = cy$mean, Rmin = cy$Rmin,
                               phi = phi, tau = tau, mu = 0), ref, label = paste(mod, "0"))
    }
})

test_that("get_rates for model k is the closed form for any mu", {
    cy <- cyc("k", mu = 0.1)
    a <- (cy$Rmax - cy$Rmin)/cy$mean
    x <- get_rates(model = "k", a = a, R = cy$mean, phi = phi, tau = tau, mu = 0.1)
    expect_equal(x$k, k, tolerance = 1e-6)
    expect_equal(x$dr, dr, tolerance = 1e-6)
})

test_that("get_rmean and get_ramp without growth rate equal the closed forms", {
    taus <- c(1, 4, 8)
    for ( mod in c("k", MODS) ) {
        expect_equal(get_rmean(k = k, k0 = k0, gamma = dr, phi = phi, tau = taus, model = mod),
                     get_rmean_nogrowth(k = k, k0 = k0, gamma = dr, phi = phi, tau = taus,
                                        model = mod), tolerance = 1e-8, label = mod)
        expect_equal(get_rmean(k = k, k0 = k0, dr = dr, mu = NA, phi = phi, tau = taus, model = mod),
                     get_rmean_nogrowth(k = k, k0 = k0, dr = dr, mu = 0, phi = phi, tau = taus,
                                        model = mod), tolerance = 1e-8, label = mod)
        for ( rel in c(TRUE, FALSE) )
            expect_equal(get_ramp(k = k, k0 = k0, gamma = dr, phi = phi, tau = taus,
                                  relative = rel, model = mod),
                         get_ramp_nogrowth(k = k, k0 = k0, gamma = dr, phi = phi, tau = taus,
                                           relative = rel, model = mod),
                         tolerance = 1e-8, label = paste(mod, rel))
    }
})

test_that("get_rates with mixed growth rates: NA counts as 0", {
    cy0 <- cyc("k_dr", mu = 0); cy1 <- cyc("k_dr", mu = 0.1)
    r <- get_rates(model = "k_dr", a = c((cy0$Rmax - cy0$Rmin)/cy0$mean, (cy1$Rmax - cy1$Rmin)/cy1$mean),
                   R = c(cy0$mean, cy1$mean), phi = phi, tau = tau, mu = c(NA, 0.1))
    expect_equal(r$dr, c(dr, dr), tolerance = 1e-6)
    expect_equal(r$k, c(k, k), tolerance = 1e-6)
})

test_that("get_rates_nogrowth and get_rates take a growth rate per gene for model k", {
    mus <- c(0.05, 0.1, NA)
    cy <- get_rcycle(k = k, dr = dr, mu = ifelse(is.na(mus), 0, mus), phi = phi, tau = tau,
                     model = "k")
    a <- (cy$Rmax - cy$Rmin)/cy$mean
    x <- get_rates_nogrowth(model = "k", a = a, R = cy$mean, phi = phi, tau = tau, mu = mus)
    expect_equal(x$dr, rep(dr, 3), tolerance = 1e-6)
    expect_equal(x$k, rep(k, 3), tolerance = 1e-6)
    y <- get_rates(model = "k", a = a, R = cy$mean, phi = phi, tau = tau, mu = mus)
    expect_equal(y, x)
})
