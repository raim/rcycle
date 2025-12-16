#!/usr/bin/Rscript --vanilla

library(rcycle)
library(unittest, quietly = TRUE)

if ( !interactive() )
    options(warn=2, error = function() {
        sink(stderr()) ; traceback(3) ; q(status = 1) })


#. temp requirements
source('/home/raim/programs/rcycle/R/models.R')
minmax <- function(x, ...) (x-min(x, ...))/(max(x, ...)-min(x, ...))
library(segmenTools)
library(deSolve)
out.path <- '/home/raim/programs/rcycle/vignettes'

## TODO:
## * move to vignettes and instead test pw_fourier, pw_simple, pw_sawtooth,

## average RP parameters (from chemostatData)
k <- 263.9
dr <- 1.7
gamma <- dr+mu
k0 <- 10

## osci and growth params
stime <- seq(0, 3, .01)
phocs <- c(.01, seq(.2, .8, .2), .99)
tosc <- 1
mu <- 0

plotdev(file.path(out.path, 'pwm_cartoon'), type='pdf', width=5, height=2)
par(mfcol=c(1,2), mai=c(.35,.5,.35,.1), mgp=c(.5,0,0))
par(mai=c(.35,.5,.35,.4))
plot(1, col=NA, xlim=range(stime), ylim=c(0, length(phocs)), axes=FALSE,
     xlab='time', ylab=NA)
text(mean(stime), length(phocs)+1,
     labels='Alternating pulse waves\nof gene expression:', font=2, xpd=TRUE)
mtext('growth rate', 2, 1)
arrows(x0=-.2, y0=.5, y1=length(phocs)-.5, xpd=TRUE, length=.1)
arrows(x0=.5, x1=max(stime)-.5, y0=-.7, xpd=TRUE, length=.1)
for ( i in seq_along(phocs) ) {
    phoc <- phocs[i]
    ploc <- 1 - phoc
    
    hocon <- as.numeric(pw_fourier(t=stime, k=1,
                                   phi=phoc, tosc=tosc, theta=0*pi,
                                   N=1e3)>.5)
    locon <- as.numeric(pw_fourier(t=stime, k=1,
                                   phi=ploc, tosc=tosc, theta=pi,
                                   N=1e3)>.5)
    
    lines(stime, .8*hocon+i-1, type='l', col=2)
    lines(stime, .8*locon+i-1, col=4, lty=3)
}
text(x=rep(max(stime), length(phocs)), y=seq_along(phocs)-.5,
     labels=phocs, col=2, xpd=TRUE, pos=4)
##text(x=max(stime), y=-.75, labels=expression(phi), col=2, xpd=TRUE, cex=1.2,
##     pos=4)
text(x=max(stime)+1, y=-1.5, labels=expression(atop(duty), cycle~phi),
     col=2, xpd=TRUE, cex=1, pos=2)


par(mai=c(.35,.5,.35,.1))
plot(phocs, col=NA, xlim=c(0, length(phocs)), ylim=c(0, 1), axes=FALSE,
     ylab=NA, xlab='growth rate')
mtext('proteome fraction', 2, 1)
arrows(x0=-.2, y0=.05, y1=.95, xpd=TRUE, length=.1)
arrows(x0=.5, x1=length(phocs)-.5, y0=-.1, xpd=TRUE, length=.1)
lines(phocs, col=2, type='b')
lines(1-phocs, col=4, type='b', lty=3)
text(length(phocs)/2, 1.15,
     labels='Continuous output:', font=2, xpd=TRUE)
dev.off()


###  CARTOON FOR PIECEWISE ODE SOLUTIONS

state <- c(R=0)
tosc <- 3
thoc <- 2
mu <- .1

## NOTE: parameters for average RP from pwmavg.R/pwmode.R
params <- c(k=263.9,
            k0=10,
            dr=1.7,
            nb=3,
            dp=0.069,
            ell=250,
            mu=mu)
times <- seq(0,50*tosc, tosc/200)
## start pulse at time 0
pulses <- pwm_simple(time=times, tosc=tosc, phi=thoc/tosc, theta = 0)
hocf <- approxfun(x = times, y = pulses,
                  method = "constant", rule = 2) # HOC pulse
outb <- ode(y = c(R=rmean(dr=params['dr'],
                          k=params['k'],
                          phi=thoc/tosc,
                          tau=tosc,
                          mu=params['mu'], model = 'k_dr')),
            times = times, func = pwmode_k_dr, 
            parms = params, hocf = hocf)
outt <- ode(y = c(R=rmean(dr=params['dr'],
                          k=params['k'],
                          phi=thoc/tosc,
                          tau=tosc,
                          mu=params['mu'],
                          model = 'k')),
            times = times, func = pwmode_k, 
            parms = params, hocf = hocf)
outd <- ode(y = c(R=rmean(dr=params['dr'],
                          k=params['k'],
                          phi=thoc/tosc,
                          tau=tosc,
                          mu=params['mu'],
                          model = 'dr')),
            times = times, func = pwmode_dr, 
            parms = params, hocf = hocf)
outk0 <- ode(y = c(R=rmean(dr=params['dr'],
                           k=params['k'],
                           k0=params['k0'],
                           phi=thoc/tosc,
                           tau=tosc,
                           mu=params['mu'],
                           model = 'k_dr_k0')),
             times = times, func = pwmode_dr, 
             parms = params, hocf = hocf)

pulse <- cbind(time=times, kappa=pulses)

## cut start for proper norm between R0 and R1
outt[outt[,1]<10*tosc,2] <- mean(outt[,2])
outb[outb[,1]<10*tosc,2] <- mean(outb[,2])
outd[outd[,1]<10*tosc,2] <- mean(outd[,2])
outk0[outk0[,1]<10*tosc,2] <- mean(outk0[,2])


plotdev(file.path(out.path, 'pwm_cartoon_piecewise'),
        type='pdf', width=3, height=2)
par(mai=c(.35,.35,.1,.35), mgp=c(1.4,0.3,0), tcl=-.25)
plot(outt[,1], minmax(outt[,2]),
     xlim=c(45*tosc-.05*tosc,
            46*tosc+.05*tosc),
     type='l', axes=FALSE,
     xlab='', ylab='')
lines(outb[,1], minmax(outb[,2]), col=2, lty=2)
lines(outd[,1], minmax(outd[,2]), col=3, lty=3)
lines(outk0[,1], minmax(outk0[,2]), col=4, lty=4)
axis(2, at=0:1, labels=expression(R[0], R[1]), las=2, col=NA)
axis(4, at=0:1, labels=expression(R[tau], R[1]), las=2, col=NA)
axis(1, at=45*tosc, labels=0, col=NA) 
axis(1, at=45*tosc+thoc, labels=expression(tau[on]), col=NA) 
axis(1, at=46*tosc, labels=expression(tau), col=NA)
abline(h=0:1, lty=3, col=8)
abline(v=seq(0:100)*tosc, lty=3, col=8)
abline(v=seq(0:100)*tosc+thoc, lty=3, col=8)
                                        #lines(pulse, col='gray')
arrows(y0=0, x0=45*tosc, x1=45*tosc+thoc,
       code=3, length=.1, col=8, lwd=1.5)
arrows(y0=0, x0=45*tosc+thoc, x1=46*tosc,
       code=3, length=.1, col=8, lwd=1.5)
if (FALSE )
    legend('topright', legend=c('const. deg.',
                                'anti-phasic deg.'), col=1:2, lty=1:2)
dev.off()
