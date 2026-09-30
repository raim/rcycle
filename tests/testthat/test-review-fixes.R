## Tests for four fixes from the code review (2026-09-30):
## 1. get_rates_nogrowth, get_times_nogrowth, get_ramp_nogrowth: the default
##    model vector failed in if(); the first model ("k") is used, as in
##    get_rmean_nogrowth and get_pmean
## 2. evaluate_order read phases$x$order from the prcomp matrix; the cohort
##    order is in phases$x.phase$order
## 3. get_cohorts, cohort2clustering: no global COHORTS/genes; n is used,
##    default the maximal index (list) or the number of columns (matrix)
## 4. center_phase wrapped only once; angles beyond -3*pi:3*pi stayed out of
##    range

library(testthat)
library(rcycle)

## 1. default model = first model, "k"
test_that("the *_nogrowth functions use the first model by default", {
    gamma <- 1; phi <- 0.3; tau <- 4; k <- 10
    R <- get_rmean_nogrowth(k = k, gamma = gamma, phi = phi, tau = tau, model = "k")
    a <- get_ramp_nogrowth(gamma = gamma, phi = phi, tau = tau, k = k, model = "k")
    expect_equal(get_ramp_nogrowth(gamma = gamma, phi = phi, tau = tau, k = k), a)
    expect_equal(get_rates_nogrowth(a = a, R = R, phi = phi, tau = tau),
                 get_rates_nogrowth(model = "k", a = a, R = R, phi = phi, tau = tau))
    expect_equal(get_rates_nogrowth(a = a, R = R, phi = phi, tau = tau)$k, k,
                 tolerance = 1e-6)
    expect_equal(get_times_nogrowth(a = a, R = R, k = k, gamma = gamma),
                 get_times_nogrowth(model = "k", a = a, R = R, k = k, gamma = gamma))
    expect_equal(unlist(get_times_nogrowth(a = a, R = R, k = k, gamma = gamma)),
                 c(phi = phi, tau = tau), tolerance = 1e-6)
})

## 2. evaluate_order
test_that("evaluate_order compares the cohort phase order with the row order", {
    ## six cohorts peaking in turn along a circle of cells
    ph <- seq(-pi, pi, length.out = 61)[-61]
    st <- t(sapply(0:5, function(i) 2 + cos(ph - i*pi/3)))
    rownames(st) <- paste0("c", 1:6)
    phases <- get_pseudophase(st, verb = 0)
    d <- evaluate_order(phases)
    expect_equal(d, state_order_distance(
        reference = rownames(phases$x),
        test = rownames(phases$x)[phases$x.phase$order]))
    expect_true(d %in% c(0, 4)) # in order, or reversed (6 letters: distance 4)
})

## 3. cohorts
test_that("get_cohorts and cohort2clustering use n, not globals", {
    coh <- list(a = 1:2, b = 4, c = 5:6)
    m <- get_cohorts(coh, n = 7)
    expect_equal(dim(m), c(3, 7))
    expect_warning(m6 <- get_cohorts(coh), "maximal index")
    expect_equal(dim(m6), c(3, 6))
    cls <- c("a", "a", "na", "b", "c", "c", "na")
    expect_equal(cohort2clustering(coh, n = 7), cls)
    expect_equal(cohort2clustering(m), cls)          # n from ncol
    expect_equal(cohort2clustering(m, n = 7), cls)
    expect_warning(cls6 <- cohort2clustering(coh), "maximal index")
    expect_equal(cls6, cls[1:6])
    expect_error(cohort2clustering(list(a = 1:2, b = 2:3), n = 3), "overlapping")
})

## 4. center_phase
test_that("center_phase brings any angle into -pi:pi", {
    x <- c(-10, -3*pi - 0.1, -pi, -1, 0, 1, pi, 3*pi + 0.1, 10, 100)
    y <- center_phase(x)
    expect_true(all(y >= -pi & y < pi))
    expect_equal(cos(y), cos(x))
    expect_equal(sin(y), sin(x), tolerance = 1e-12)
    ## within -3*pi:3*pi, unchanged from the previous two-step version
    z <- seq(-3*pi, 3*pi - 1e-9, length.out = 1001)
    old <- z
    old[old <  pi] <- old[old <  pi] + 2*pi
    old[old >= pi] <- old[old >= pi] - 2*pi
    expect_identical(center_phase(z), old)
})
