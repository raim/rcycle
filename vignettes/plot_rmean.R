#!/usr/bin/Rscript --vanilla

library(rcycle)

#. temp requirements
library(segmenTools)
source('/home/raim/programs/rcycle/R/models.R')
out.path <- '/home/raim/programs/rcycle/vignettes'

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

## TODO: test the hyperbolic term
if ( FALSE ) {

    y <- gamma*tau*(1-phis)
    ey <- pracma::coth(y/2)

    plot(phis, ey)
    x <- phis^2*k*tau/2
    b <- phis*k/gamma
    plot(phis, b, ylim=c(0,500), type='l')
    lines(phis, x)
    lines(phis, x*ey)
}


## relative amplitudes dependence on period
## NOTE: higher abundance with increasing period!
models <- c('k', 'dr', 'k_dr', 'k_dr_k0')
tmns <- matrix(NA, nrow=length(taus), ncol=length(models))
colnames(tmns) <- models
for ( mod in models ) {
    tmns[,mod] <- get_rmean(k=k, gamma=gamma, k0=k0, phi=phi, tau=taus,
                            model = mod, use.coth = TRUE)
}

plotdev(file.path(out.path, 'pwm_rmean_tau'),
        type='pdf', width=3, height=3)
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(taus, tmns, type='l', lty=1, col=1:ncol(tmns),
        ylim=c(0, max(tmns, na.rm=TRUE)),        
        xlab=expression(period~tau), ylab=axis_labels['rmean'])
legend('topleft', colnames(pmns), col=1:ncol(pmns), lty=1, bty='n')         
dev.off()

## relative amplitudes dependence on duty cycle phi
## NOTE: average abundance explodes towards phi->1.
pmns <- matrix(NA, nrow=length(phis), ncol=length(models))
colnames(pmns) <- models
for ( mod in models ) {

    pmns[,mod] <- get_rmean(k=k, gamma=gamma, k0=k0, phi=phis, tau=tau,
                            model = mod, use.coth = TRUE)
}

plotdev(file.path(out.path, 'pwm_rmean_phi'),
        type='pdf', width=3, height=3)
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(phis, pmns, type='l', lty=1, col=1:ncol(pmns),
        ylim=c(0, 3*max(pmns[,'k'], na.rm=TRUE)),
        xlab=expression(duty~cycle~phi), ylab=axis_labels['rmean'])
legend('topleft', colnames(pmns), col=1:ncol(pmns), lty=1, bty='n')         
dev.off()

plotdev(file.path(out.path, 'pwm_rmean_phi_log'),
        type='pdf', width=3, height=3)
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(phis, log10(pmns), type='l', lty=1, col=1:ncol(pmns),
        ylim=c(0, log10(max(pmns[is.finite(pmns)], na.rm=TRUE))),
        axes=FALSE, xlab=expression(duty~cycle~phi), ylab=axis_labels['rmean'])
axis(1)
segmenTools::logaxis(2)
legend('bottomright', colnames(pmns), col=1:ncol(pmns), lty=1, bty='n')         
dev.off()


