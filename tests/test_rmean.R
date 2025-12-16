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

## test the hyperbolic term
y <- gamma*tau*(1-phis)
ey <- pracma::coth(y/2)

plot(phis, ey)

x <- phis^2*k*tau/2

b <- phis*k/gamma
plot(phis, b, ylim=c(0,500), type='l')
lines(phis, x)
lines(phis, x*ey)


models <- c('k', 'dr', 'k_dr', 'k_dr_k0')
tmns <- matrix(NA, nrow=length(taus), ncol=length(models))
colnames(tmns) <- models
for ( mod in models ) {
    tmns[,mod] <- rmean(k=k, gamma=gamma, k0=k0, phi=phi, tau=taus,
                        model = mod, use.coth = TRUE)
}

matplot(taus, tmns, type='l', lty=1, col=1:ncol(tmns),
        ylim=c(0, max(tmns, na.rm=TRUE)),        
        xlab=expression(period~tau), ylab=axis_labels['rmean'])
legend('topleft', colnames(pmns), col=1:ncol(pmns), lty=1)         

pmns <- matrix(NA, nrow=length(phis), ncol=length(models))
colnames(pmns) <- models
for ( mod in models ) {

    pmns[,mod] <- rmean(k=k, gamma=gamma, k0=k0, phi=phis, tau=tau,
                        model = mod, use.coth = TRUE)
}

matplot(phis, pmns, type='l', lty=1, col=1:ncol(pmns),
        ylim=c(0, max(pmns[is.finite(pmns)], na.rm=TRUE)),
        xlab=expression(duty~cycle~phi), ylab=axis_labels['rmean'])
matplot(phis, log10(pmns), type='l', lty=1, col=1:ncol(pmns),
        ylim=c(0, log10(max(pmns[is.finite(pmns)], na.rm=TRUE))),
        axes=FALSE, xlab=expression(duty~cycle~phi), ylab=axis_labels['rmean'])
axis(1)
segmenTools::logaxis(2)
legend('bottomright', colnames(pmns), col=1:ncol(pmns), lty=1)         



## Test equivalence of expm1 and coth-based implementations
## coth(y/2) = (e^y+1)/(e^y-1)
ok(ut_cmp_equal(rmean(k=k, dr=dr, mu=mu, phi=.5, tau=taus,
                      model = 'dr', use.coth = TRUE),
                rmean(k=k, dr=dr, mu=mu, phi=.5, tau=taus,
                      model = 'dr', use.coth = FALSE)), "two numbers")
