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
k0 <- 10

## osci and growth params
phi <- .5
tau <- 2

## vary over tau and phi
taus <- seq(0,8, .1)
phis <- 0:100/100

## relative amplitudes dependence on period tau
models <- c('k', 'dr', 'k_dr', 'k_dr_k0')
tmns <- matrix(NA, nrow=length(taus), ncol=length(models))
colnames(tmns) <- models
for ( mod in models ) {
    tmns[,mod] <- get_ramp(gamma=gamma, phi=phi, tau=taus,
                           k=k, k0=k0, # only required with basal expression!
                           model = mod, relative=TRUE)
}

plotdev(file.path(out.path, 'pwm_ramp_tau'),
        type='pdf', width=3, height=3)
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(taus, tmns, type='l', lty=1, col=1:ncol(tmns),
        ylim=c(0,10), #range(tmns, na.rm=TRUE),        
        xlab=expression(period~tau), ylab=axis_labels['rampr'])
legend('topleft', colnames(pmns), col=1:ncol(pmns), lty=1)         
dev.off()

## relative amplitudes dependence on duty cycle phi
pmns <- matrix(NA, nrow=length(phis), ncol=length(models))
colnames(pmns) <- models
for ( mod in models ) {
    pmns[,mod] <- get_ramp(gamma=gamma, phi=phis, tau=tau,
                           k=k, k0=k0, # only required with basal expression!
                           model = mod, relative=TRUE)
}

plotdev(file.path(out.path, 'pwm_ramp_phi'),
        type='pdf', width=3, height=3)
par(mai=c(.5,.5,.15,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(phis, pmns, type='l', lty=1, col=1:ncol(pmns),
        ylim=c(0, max(pmns[is.finite(pmns)], na.rm=TRUE)),
        xlab=expression(duty~cycle~phi), ylab=axis_labels['rampr'])
matplot(phis, log10(pmns), type='l', lty=1, col=1:ncol(pmns),
        ylim=c(log10(0.03), log10(max(pmns[is.finite(pmns)], na.rm=TRUE))),
        axes=FALSE, xlab=expression(duty~cycle~phi), ylab=axis_labels['rampr'])
axis(1)
segmenTools::logaxis(2)
legend('bottomleft', colnames(pmns), col=1:ncol(pmns), lty=1, bty='n')         
dev.off()

