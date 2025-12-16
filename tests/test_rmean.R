#!/usr/bin/Rscript --vanilla

library(rcycle)
library(unittest, quietly = TRUE)

if ( !interactive() )
    options(warn=2, error = function() {
        sink(stderr()) ; traceback(3) ; q(status = 1) })

#. temp requirements
library(segmenTools)
source('/home/raim/programs/rcycle/R/models.R')
out.path <- '/home/raim/programs/rcycle/vignettes'

## TODO:
## * add and load common parameter set,
## * Generate more tests of model results:
##   This currently only tests coth vs. expm1 implementations,
##   which is trivial,

## average RP parameters (from chemostatData)
k <- 263.9
dr <- 1.7
mu <- 0
gamma <- dr+mu
k0 <- k/10

## osci and growth params
phi <- .5
tau <- 2

## vary over tau and phi
taus <- seq(0,8, .1)
phis <- 0:100/100



## Test equivalence of expm1 and coth-based implementations
## coth(y/2) = (e^y+1)/(e^y-1)

models <- c('k', 'dr', 'k_dr', 'k_dr_k0')
for ( mod in models ) {
    ok(ut_cmp_equal(get_rmean(k=k, k0=k0, dr=dr, mu=mu, phi=.5, tau=taus,
                              model = mod, use.coth = TRUE),
                    get_rmean(k=k, k0=k0, dr=dr, mu=mu, phi=.5, tau=taus,
                              model = mod, use.coth = FALSE)), "two numbers")
}
