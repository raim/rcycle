#!/usr/bin/Rscript --vanilla

library(rcycle)
library(testthat)

#. temp requirements
library(segmenTools)
source('/home/raim/programs/rcycle/R/models.R')
out.path <- '/home/raim/programs/rcycle/vignettes'


## TODO:
## * generate more tests of model results.
##   currently tests consistency of amplitude and mean calc.

## average RP parameters (from chemostatData)
k <- 263.9
dr <- 1.7
mu <- 0
gamma <- dr+mu
k0 <- 10

## osci and growth params
phi <- .5
tau <- 2

## vary over tau and phi
taus <- seq(0,8, .1)
phis <- 0:100/100

## relative amplitudes dependence on period tau
models <- c('k', 'dr', 'k_dr', 'k_dr_k0')

## Test equivalence of relative amplitude rampr == ramp/rmean
for ( mod in models ) {

    test_that("relative and absolute amplitudes agree", {
        ampr <- get_ramp(gamma=gamma, phi=phis,
                         tau=tau, model = mod, relative=TRUE,
                         k=k, k0=k0)
        amp <- get_ramp(gamma=gamma, phi=phis,
                        tau=tau, model = mod, relative=FALSE,
                        k=k, k0=k0)
        mn <- get_rmean(k=k, k0=k0, gamma=gamma, phi=phis, tau=tau, model = mod)
    
        ## should be equiv.
        ##plot(ampr, amp/mn); abline(a=0, b=1)

        ampm <- (amp/mn)
        nan <- !is.na(ampm)
        
        ## TODO: test expected number of NA
        ##sum(!nan)
        expect_equal(ampr[nan], ampm[nan], tolerance = 1e-10)
    })
}
