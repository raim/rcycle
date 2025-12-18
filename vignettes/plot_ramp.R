#!/usr/bin/Rscript --vanilla

library(rcycle)

#. temp requirements
library(segmenTools)
source('/home/raim/programs/rcycle/R/models.R')
out.path <- '/home/raim/programs/rcycle/vignettes'


## TODO:
## * move plots to vignette,
## * generate more tests of model results,
## * add transcriptome osci conditions!

W <- H <- 2.5

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
taus <- seq(0,7.5, .1)
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
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(taus, tmns, type='l', lty=1, col=1:ncol(tmns),
        ylim=c(0,6), #range(tmns, na.rm=TRUE),        
        xlab=expression(period~tau/h), ylab=axis_labels['rampr'])
legend('topleft', colnames(tmns), col=1:ncol(tmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)
axis(4, labels=FALSE)
mtext(bquote(phi==.(phi)), 3, 0)
dev.off()

plotdev(file.path(out.path, 'pwm_ramp_tau_log'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(taus, log10(tmns), type='l', lty=1, col=1:ncol(tmns),
        axes=FALSE,
        xlab=expression(period~tau/h), ylab=axis_labels['rampr'])
logaxis(2)
axis(1)
logaxis(4, labels=FALSE)
box()
legend('topright', colnames(tmns), col=1:ncol(tmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
mtext(bquote(phi==.(phi)), 3, 0)
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
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
ppmns <- pmns
ppmns[ppmns<0] <- NA
matplot(phis, ppmns, type='l', lty=1, col=1:ncol(pmns),
        ylim=c(0, 6), #max(pmns[is.finite(pmns)], na.rm=TRUE)),
        xlab=expression(duty~cycle~phi), ylab=axis_labels['rampr'])
legend('topright', colnames(pmns), col=1:ncol(pmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)
axis(4, labels=FALSE)
mtext(bquote(tau==.(tau)~h), 3, 0)
dev.off()

plotdev(file.path(out.path, 'pwm_ramp_phi_log'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.5,.25,.15), mgp=c(1.4,0.3,0), tcl=-.25)
matplot(phis, log10(pmns), type='l', lty=1, col=1:ncol(pmns),
        ylim=c(log10(0.03), log10(max(pmns[is.finite(pmns)], na.rm=TRUE))),
        axes=FALSE, xlab=expression(duty~cycle~phi), ylab=axis_labels['rampr'])
axis(1)
logaxis(2)
logaxis(4, labels=FALSE)
box()
legend('bottomleft', colnames(pmns), col=1:ncol(pmns), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
mtext(bquote(tau==.(tau)~h), 3, 0)
dev.off()

