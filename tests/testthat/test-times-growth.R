## Exact duty cycle and period with dilution in both phases (get_times): the
## inverse of get_rcycle for phi and tau, given the rates; without growth and
## for model k the closed forms (get_times_nogrowth).

library(testthat)
library(rcycle)

## k, k0, dr, mu, phi, tau; the last set is in the reversed regime for k_dr_k0
PARS <- list(c(k = 10, k0 = 2, dr = 1, mu = 0.1, phi = 0.3, tau = 4),
             c(k = 10, k0 = 5, dr = 3, mu = 0.2, phi = 0.6, tau = 1),
             c(k = 5, k0 = 20, dr = 0.3, mu = 0.05, phi = 0.2, tau = 8),
             c(k = 1, k0 = 20, dr = 0.3, mu = 0.2, phi = 0.3, tau = 6))
fwd <- function(p, mod, mu = p[["mu"]])
    get_rcycle(k = p[["k"]], k0 = p[["k0"]], dr = p[["dr"]], mu = mu,
               phi = p[["phi"]], tau = p[["tau"]], model = mod)

test_that("get_times inverts get_rcycle for all models, incl. the reversed regime", {
    for ( i in seq_along(PARS) ) for ( mod in c("k", "dr", "k_dr", "k_dr_k0") ) {
        p <- PARS[[i]]; cy <- fwd(p, mod)
        tm <- get_times(model = mod, a = (cy$Rmax - cy$Rmin)/cy$mean, R = cy$mean,
                        k = p[["k"]], k0 = p[["k0"]], dr = p[["dr"]], mu = p[["mu"]])
        lab <- paste(mod, paste(names(p), p, sep = "=", collapse = " "))
        expect_equal(tm$phi, p[["phi"]], tolerance = 1e-6, label = paste(lab, "phi"))
        expect_equal(tm$tau, p[["tau"]], tolerance = 1e-6, label = paste(lab, "tau"))
    }
})

test_that("get_times for k_dr_k0 uses Rmin when k0 is not given", {
    for ( p in PARS ) {
        cy <- fwd(p, "k_dr_k0")
        tm <- get_times(model = "k_dr_k0", a = (cy$Rmax - cy$Rmin)/cy$mean, R = cy$mean,
                        Rmin = cy$Rmin, k = p[["k"]], dr = p[["dr"]], mu = p[["mu"]])
        expect_equal(tm$phi, p[["phi"]], tolerance = 1e-6)
        expect_equal(tm$tau, p[["tau"]], tolerance = 1e-6)
    }
})

test_that("get_times without growth rate uses the closed forms", {
    p <- PARS[[1]]
    for ( mod in c("k", "dr", "k_dr", "k_dr_k0") ) {
        cy <- fwd(p, mod, mu = 0)
        a <- (cy$Rmax - cy$Rmin)/cy$mean
        ref <- get_times_nogrowth(model = mod, a = a, R = cy$mean, Rmin = cy$Rmin,
                                  k = p[["k"]], gamma = p[["dr"]])
        for ( mu in list(NULL, NA, 0) ) {
            args <- list(model = mod, a = a, R = cy$mean, Rmin = cy$Rmin, k = p[["k"]],
                         dr = p[["dr"]])
            if ( !is.null(mu) ) args$mu <- mu
            expect_equal(do.call(get_times, args), ref, label = paste(mod, format(mu)))
        }
        ## the closed form recovers the times at mu = 0
        expect_equal(ref$phi, p[["phi"]], tolerance = 1e-6, label = mod)
        expect_equal(ref$tau, p[["tau"]], tolerance = 1e-5, label = mod)
    }
})

test_that("get_times_nogrowth with mu = NA does not return NA", {
    p <- PARS[[1]]; cy <- fwd(p, "k_dr", mu = 0)
    x <- get_times_nogrowth(model = "k_dr", a = (cy$Rmax - cy$Rmin)/cy$mean, R = cy$mean,
                            k = p[["k"]], dr = p[["dr"]], mu = NA)
    expect_equal(x$tau, p[["tau"]], tolerance = 1e-5)
})

test_that("get_times is vectorised, takes gamma for dr, and returns NA without a solution", {
    p <- PARS[[1]]
    cy <- get_rcycle(k = p[["k"]], dr = p[["dr"]], mu = p[["mu"]], phi = c(0.2, 0.5),
                     tau = c(2, 6), model = "k_dr")
    tm <- get_times(model = "k_dr", a = (cy$Rmax - cy$Rmin)/cy$mean, R = cy$mean,
                    k = p[["k"]], gamma = p[["dr"]] + p[["mu"]], mu = p[["mu"]])
    expect_equal(tm$phi, c(0.2, 0.5), tolerance = 1e-6)
    expect_equal(tm$tau, c(2, 6), tolerance = 1e-6)
    ## a mean no duty cycle can produce
    tm <- get_times(model = "k_dr", a = 0.5, R = 1e6, k = p[["k"]], dr = p[["dr"]], mu = p[["mu"]])
    expect_true(is.na(tm$tau))
})
