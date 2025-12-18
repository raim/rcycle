#!/usr/bin/Rscript --vanilla

library(rcycle)
library(testthat)


#. temp requirements
source('/home/raim/programs/rcycle/R/models.R')
minmax <- function(x, ...) (x-min(x, ...))/(max(x, ...)-min(x, ...))
library(segmenTools)
library(deSolve)
out.path <- '/home/raim/programs/rcycle/vignettes'

## TODO:
## * Write proper tests for pulse waves, 
##   approx. function, and ODE models, numeric vs. analytic, 
## * Decide on which pulse wave functions to keep, or analyze.

## average RP parameters (from chemostatData)
k <- 263.9
dr <- 1.7
gamma <- dr+mu
k0 <- 10

## osci and growth params
stime <- seq(0, 3, .01)
phocs <- c(.01, seq(.2, .8, .2), .99)
tosc <- 2
thoc <- .5
mu <- 0


## test pulse wave functions
pwp <- pwm_simple(time=stime, tosc=tosc, phi=thoc/tosc, theta = 0)
pwf <- pw_fourier(t=stime, k=1,
                  phi=thoc/tosc, tosc=tosc, theta=0*pi,
                  N=1e3)
pws <- pw_sawtooth(t=stime, k=1,
                   thoc=thoc, tosc=tosc, alpha=0)

## * test analytic vs. ODE.
if ( interactive() ) {
    plot(stime, pwp, type='l')
    lines(stime, pwf, col=2)
    lines(stime, pws, col=3)
}

## ODE
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


