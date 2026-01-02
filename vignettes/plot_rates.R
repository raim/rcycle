#!/usr/bin/Rscript --vanilla

library(rcycle)

#. temp requirements
library(segmenTools)
source('/home/raim/programs/rcycle/R/models.R')
out.path <- '/home/raim/programs/rcycle/vignettes'

## TODO:
## * WHY DO WE GET A U SHAPE FOR k(phi) for the constant model?
## * implement all and compare coth and expm1 versions.
## * find better regions for k_dr_k0!

W <- H <- 2.5

## average RP parameters 
rmax <- 120
rmin <- 30
rmean <- mean(c(rmax, rmin))
ramp <- rmax - rmin

## osci and growth params
phi <- .5
tau <- 2
mu <- 0

## vary over tau and phi
taus <- seq(0,7.5, .1)
phis <- 1:100/100


models <- c('k', 'dr', 'k_dr', 'k_dr_coth', 'k_dr_k0', 'k_dr_k0_coth') 

mrates <- list()
krng <- drrng <- c()
for ( mod in models ) {

    a <- seq(.5,5,.4)
    
    mrates[[mod]] <-
        get_rates(model = mod,
                  A=ramp, R=rmean, Rmin=rmin, 
                  phi=phi, tau=tau, mu=mu,
                  lower = 1e-6, upper = 1e4, verb = 0)
    
    tmp <- mrates[[mod]][,'k']
    krng <- range(c(krng, tmp[is.finite(tmp)]), na.rm=TRUE)
    tmp <- mrates[[mod]][,'dr']
    drrng <- range(c(drrng, tmp[is.finite(tmp)]), na.rm=TRUE)

}


prates <- list()
krng <- drrng <- c()
for ( mod in models ) {
    
    prates[[mod]] <-
        get_rates(model = mod,
                  A=ramp, R=rmean, Rmin=rmin, Rmax=rmax,
                  phi=phis, tau=tau, mu=mu,
                  lower = 1e-6, upper = 1e4, verb = 0)
    
    tmp <- prates[[mod]][,'k']
    krng <- range(c(krng, tmp[is.finite(tmp)]), na.rm=TRUE)
    tmp <- prates[[mod]][,'dr']
    drrng <- range(c(drrng, tmp[is.finite(tmp)]), na.rm=TRUE)

}

plotdev(file.path(out.path, 'pwm_rates_k_phi'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.55,.25,.15), mgp=c(1.3,0.3,0), tcl=-.25)
plot(phis, phis, ylim=log10(krng), col=NA, axes=FALSE,
     xlab=axis_labels['phi'], ylab=axis_labels['k'])
axis(1)
logaxis(2)
logaxis(4, labels=FALSE)
for ( i in seq_along(models) ) 
    lines(phis, log10(prates[[models[[i]]]][,'k']), col=i, type='l', lty=i)
mtext(bquote(tau==.(tau)~h), 3, 0)
legend('topright', models, col=1:length(models), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
dev.off()

plotdev(file.path(out.path, 'pwm_rates_dr_phi'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.55,.25,.15), mgp=c(1.3,0.3,0), tcl=-.25)
plot(phis, phis, ylim=log10(drrng), col=NA, axes=FALSE,
     xlab=axis_labels['phi'], ylab=axis_labels['gamma'])
axis(1)
logaxis(2)
logaxis(4, labels=FALSE)
for ( i in seq_along(models) )
    lines(phis, log10(prates[[models[[i]]]][,'dr']), col=i, type='l', lty=i)
mtext(bquote(tau==.(tau)~h), 3, 0)
legend('bottomright', models, col=1:length(models), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
dev.off()


trates <- list()
for ( mod in models ) {
    
    trates[[mod]] <-
        get_rates(model = mod,
                  A=ramp, R=rmean, Rmin=rmin, Rmax=rmax,
                  phi=phi, tau=taus, mu=mu,
                  lower = 1e-6, upper = 1e4, verb = 0)

    tmp <- trates[[mod]][,'k']
    krng <- range(c(krng, tmp[is.finite(tmp)]), na.rm=TRUE)
    tmp <- trates[[mod]][,'dr']
    drrng <- range(c(drrng, tmp[is.finite(tmp)]), na.rm=TRUE)
}

plotdev(file.path(out.path, 'pwm_rates_k_tau'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.55,.25,.15), mgp=c(1.3,0.3,0), tcl=-.25)
plot(taus, taus, ylim=log10(krng), col=NA, axes = FALSE,
     xlab=axis_labels['tau'], ylab=axis_labels['k'])
for ( i in seq_along(models) ) 
    lines(taus, log10(trates[[models[[i]]]][,'k']), col=i, type='l', lty=i)
axis(1)
logaxis(2)
logaxis(4, labels=FALSE)
mtext(bquote(phi==.(phi)), 3, 0)
legend('topright', models, col=1:length(models), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
dev.off()

plotdev(file.path(out.path, 'pwm_rates_dr_tau'),
        type='pdf', width=W, height=H)
par(mai=c(.5,.55,.25,.15), mgp=c(1.3,0.3,0), tcl=-.25)
plot(taus, taus, ylim=log10(drrng), col=NA, axes=FALSE,
     xlab=axis_labels['tau'], ylab=axis_labels['gamma'])
for ( i in seq_along(models) ) 
    lines(taus, log10(trates[[models[[i]]]][,'dr']), col=i, type='l', lty=i)
axis(1)
logaxis(2)
logaxis(4, labels=FALSE)
mtext(bquote(phi==.(phi)), 3, 0)
legend('bottomleft', models, col=1:length(models), lty=1, bty='n',
       seg.len=.5, y.intersp=.75)         
dev.off()
