#!/usr/bin/Rscript --vanilla

library(rcycle)
library(unittest, quietly = TRUE)

if ( !interactive() )
    options(warn=2, error = function() {
        sink(stderr()) ; traceback(3) ; q(status = 1) })

source('/home/raim/programs/rcycle/R/models.R')

## TODO:
## * move plots to vignette,
## * generate more tests of model results.

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

pamp <- ramp(gamma=gamma, phi=phis, tau=tau, model = 'k', relative=TRUE)

plot(phis, pamp)

tamp <- ramp(gamma=gamma, phi=phi, tau=taus, 
             k=k, k0=k0, # only required for model with basal expression!
             model = 'k_dr_k0', relative=TRUE)

plot(taus, tamp)

models <- c('k', 'dr', 'k_dr', 'k_dr_k0')
tmns <- matrix(NA, nrow=length(taus), ncol=length(models))
colnames(tmns) <- models
for ( mod in models ) {
    tmns[,mod] <- ramp(gamma=gamma, phi=phi, tau=taus,
                       k=k, k0=k0, # only required with basal expression!
                       model = mod, relative=TRUE)
}

matplot(taus, tmns, type='l', lty=1, col=1:ncol(tmns),
        ylim=c(0,10), #range(tmns, na.rm=TRUE),        
        xlab=expression(period~tau), ylab=axis_labels['rampr'])
legend('topleft', colnames(pmns), col=1:ncol(pmns), lty=1)         

pmns <- matrix(NA, nrow=length(phis), ncol=length(models))
colnames(pmns) <- models
for ( mod in models ) {
    pmns[,mod] <- ramp(gamma=gamma, phi=phis, tau=tau,
                       k=k, k0=k0, # only required with basal expression!
                       model = mod, relative=TRUE)
}

matplot(phis, pmns, type='l', lty=1, col=1:ncol(pmns),
        ylim=c(0, max(pmns[is.finite(pmns)], na.rm=TRUE)),
        xlab=expression(duty~cycle~phi), ylab=axis_labels['rampr'])
matplot(phis, log10(pmns), type='l', lty=1, col=1:ncol(pmns),
        ylim=c(log10(0.01), log10(max(pmns[is.finite(pmns)], na.rm=TRUE))),
        axes=FALSE, xlab=expression(duty~cycle~phi), ylab=axis_labels['rampr'])
axis(1)
segmenTools::logaxis(2)
legend('topright', colnames(pmns), col=1:ncol(pmns), lty=1)         


## Test equivalence of relative amplitude rampr == ramp/rmean
for ( mod in models ) {

    ampr <- ramp(gamma=gamma, phi=phis, tau=tau, model = mod, relative=TRUE,
                 k=k, k0=k0)
    amp <- ramp(gamma=gamma, phi=phis, tau=tau, model = mod, relative=FALSE,
                k=k, k0=k0)
    mn <- rmean(k=k, k0=k0, gamma=gamma, phi=phis, tau=tau, model = mod)

    ## should be equiv.
    ##plot(ampr, amp/mn); abline(a=0, b=1)

    if ( model!='k' )
        ok(ut_cmp_equal(tail(mn, 1), Inf, "two numbers"))

    nan <- is.finite(mn)
    ok(ut_cmp_equal(ampr[nan], (amp/mn)[nan], "two numbers"))
}
