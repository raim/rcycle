
## MODELS for PERIODIC GENE EXPRESSION

## ODE Models

#' ODE of pulse width-modulated transcription.
#' @export
pwmode_k <- function(time, state, parameters, hocf){
    kappa <- hocf(time)
    with(as.list(c(state, parameters)), {
        ## NOTE: compared to older functions, here
        ## kappa is 0 or 1!
        dR = kappa*k - (dr+mu)*R
        return(list(c(dR)))
    })
}

#' ODE of pulse width-modulated degradation.
#' @export
pwmode_dr <- function(time, state, parameters, hocf){
    kappa <- hocf(time)
    with(as.list(c(state, parameters)), {
        dR = k - ((1-kappa)*dr+mu)*R
        return(list(c(dR)))
    })
}

#' ODE of pulse width-modulated transcription and anti-phasic
#' degradation.
#' @export
pwmode_k_dr <- function(time, state, parameters, hocf){
    kappa <- hocf(time)
    with(as.list(c(state, parameters)), {
        dR = kappa*k - ((1-kappa)*dr+mu)*R
        return(list(c(dR)))
    })
}

#' ODE of pulse width-modulated transcription and anti-phasic degradation
#' and basal transcription.
#' @export
pwmode_k_dr_k0 <- function(time, state, parameters, hocf){
    kappa <- hocf(time)
    with(as.list(c(state, parameters)), {
        dR = kappa*k + (1-kappa)*k0 - ((1-kappa)*dr+mu)*R
        return(list(c(dR)))
    })
}


### ANALYTIC

#' Calculate mean abundance from rates and times.
#'@export
rmean <- function(k, gamma, k0, dr, mu, phi, tau,
                  model = c('k', 'dr', 'k_dr', 'k_dr_k0'),
                  use.coth = TRUE) {

    if ( missing(gamma) )
        gamma <- dr+mu

    ## base model: phi*k/gamma
    rmn <- phi*k/gamma

    ## expand, if the model is requested for a vector of periods
    if ( !missing(tau) ) 
        rmn <- rep(rmn, length(tau))

    ## add terms for all other models
    ## phi^2*k*tau/2 * coth(gamma*tau*(1-phi)/2)
    if ( model %in% c('dr', 'k_dr', 'k_dr_k0') ) {

        ## exponent
        y <- gamma*tau*(1-phi)
        if ( use.coth ) 
            term2 <- pracma::coth(y/2)
        else {
            ye <- expm1(y) # exp(y)-1
            term2 <- (ye+2)/ye
        }

        rmn <- rmn + (phi^2*k*tau/2)*term2
    }
    if ( model %in% c('dr') )
        rmn <- rmn + k/gamma
    
    if ( model %in% c('k_dr_k0') )
        rmn <- rmn + k0/gamma

    rmn    
    
}


ramp <- function(gamma, dr, mu, phi, tau, relative = TRUE, k, k0,
                 model = c('k', 'dr', 'k_dr', 'k_dr_k0'), ...) {

    if ( model %in% c('dr', 'k_dr', 'k_dr_k0') ) {

        ## calculate relative amplitude
        if ( missing(gamma) )
            gamma <- dr+mu

        gt <- gamma*tau

        term2 <- phi/2 * coth(gt*(1-phi)/2)

        term1 <- switch (model,
                         k_dr = 1/gt,
                         dr = (1-1/phi)/gt,
                         k_dr_k0 = (1-k0/(k*phi))/gt)

        ramp <- 1/(term1 + term2)

        if ( missing(k0) ) k0 <- NA

        if ( !relative ) {
            rmean <- rmean(
                k = k,
                gamma = gamma,
                phi = phi,
                tau = tau,
                k0 = k0,
                model = model, ...)
            ramp <- ramp*rmean 
        }
    } else if ( model %in% c('k') ) {

        if ( missing(gamma) )
            gamma <- dr+mu

        gt <- gamma*tau
        term1 <- 1-exp(-gt*phi)
        term2 <- 1-exp(gt*(phi-1))
        term3 <- 1-exp(-gt) 

        term <- term1*term2/term3

        if ( relative ) ramp <- term/phi
        else ramp <- term*k/gamma
    }

    ramp
    
}



get_rates <- function(model = c('k', 'dr', 'k_dr', 'k_dr_k0'),
                      ramp, rmean, rmin, rmax,
                      phi=phis, tau=tosc, mu=mu,
                      lower = 1e-6, upper = 20, ...) {

    

}

get_times <- function(model = c('k', 'dr', 'k_dr', 'k_dr_k0'),
                      k, dr, k0, lower = 1e-6, upper = 20, ...) {
}


## transcription with const. degradation,
## where x = gamma * tosc
root_k <- function(x) {
    lhs <- a * phi
    ##term1 <- (1 - exp(-phi * x)) / (1 - exp(-x))
    ##term2 <- (1 - exp(x * (phi - 1)))
    term1 <- (- expm1(-phi * x)) / (- expm1(-x))
    term2 <- ( - expm1(x * (phi - 1)))
    rhs <- term1 * term2
    return(lhs - rhs)
}

## anti-phasic transcription/degradation,
## where x = gamma * tosc
root_k_dr <- function(x) {

    lhs <- 1/a
    ##term1 <- phi/(exp(x * (1-phi))-1)
    term1 <- phi/expm1(x * (1-phi)) # = exp() + 1, where 1 cancels below
        term2 <- phi/2
    term3 <- 1/x
    rhs <- term1 + term2 + term3
    return(lhs - rhs)
}

## anti-phasic transcription/degradation and basal transcription
## where x = gamma * tosc
root_k_dr_k0 <- function(x) {

    ## requires mean(R), min(R), and either amplitde a or
    ## transcription rate k
    
    lhs <- R - Rmin

    betam1 <- expm1(x*(1-phi)) # exp() +1
    
    term1 <- (phi-1) /betam1 
    term2 <- phi/2
    term3 <- 1/x 
    
    rhs <- phi*k*tosc*(term1 + term2 + term3)

    return(lhs - rhs)
}


#' PWM of transcription, analytic solution.
#'
#' Analytic solution of the pulse-wave ODE for transcript abundance,
#' optionally extended to proteins via numeric integration.
#'
#' @param t Numeric vector of time points.
#' @param R0 Initial RNA abundance at t=0.
#' @param k Transcription rate during pulse (ON state).
#' @param gamma total turnover rate, gamma=dr+mu.
#' @param dr RNA degradation rate, required if gamma is missing.
#' @param mu Growth/dilution rate, required if gamma is missing.
#' @param phi Duty cycle (fraction of period ON, between 0 and 1).
#' @param tosc Oscillation period.
#' @param k0 Basal transcription rate (default 0).
#' @param alpha Fourier damping factor (default 0).
#' @param theta Phase shift in radians (default 0).
#' @param N Number of Fourier terms to use (default 500).
#' @param P0 protein abundance at time t[1] (!).
#' @param dp protein degradation rate.
#' @param ell transcript elongation rate.
#' @param rho translating ribosomes per mRNA.
#'
#' @return Data frame with columns:
#'   - `time`: original time points
#'   - `R`: RNA abundance
#'   - `pulse`: binary 0/1 pulse state
#'   - `P`: protein abundance (if protein params given)
#'@export
pwm_k <- function(t, R0, k, gamma, dr,  mu, k0=0, 
                  phi, tosc,
                  alpha=0, shift=0, theta=0, N=5e2,
                  P0, dp, ell, rho) {


    otime <- t # store original time before shifting
    
    ## phase shift by shifting time
    ## TODO: fix this, this simple solution shifts R0 to t=shift
    if ( shift != 0 ) {
        t <- t - shift
        warning("phase shift also shifts R0")
    }
    
    ## total turnover
    if ( missing(gamma) )
        gamma <- mu+dr

    ## angular frequency
    omega <- 2*pi/tosc

    ## exponential decay term
    EGT <- exp(-gamma*t)
    
    ## contribution from initial concentration
    r0 <- EGT*R0

    ## contribution from constant term
    rc <- (k*phi + k0)/gamma * (1-EGT)

    ## contribution from periodic term

    ## phase shift theta present?
    shifted <- !missing(theta)
    if ( shifted ) shifted <- shifted & theta!=0
    
    sm <- 0
    if ( !shifted ) 
        for ( n in 1:N ) { # NO PHASE SHIFT
            
            ra <- (1/n) * sin(n*pi*phi) * exp(-alpha*n)
            nom <- t*n*omega 
            ra <- ra * (cos(nom) + n*omega*sin(nom)/gamma - EGT)
            ra <- ra / (gamma + n^2*omega^2/gamma)
            
            sm <- sm + ra
        }
    else
        for ( n in 1:N ) { # WITH PHASE SHIFT
            
            ra <- (1/n) * sin(n*pi*phi) * exp(-alpha*n)
            tnom <- n*(t*omega - theta) 
            ra <- ra * (gamma*cos(tnom) + n*omega*sin(tnom) +
                        EGT *(-gamma*cos(n*theta) + n*omega*sin(n*theta)))
            ra <- ra / (gamma^2 + n^2*omega^2)
            
            sm <- sm + ra
        }
    rp <- 2*k/pi *sm
        
    rt <- data.frame(time=otime, R=r0+rc+rp)

    ## protein if parameters are present!
    ## TODO:
    ## * test and fix this, esp. for times >> 0!!??
    ## * for t>0 do we need to still calculate all R from t==0 and integrate?
    ## * get proper analytic solution instead of just integrating over RNA?
    if ( !missing(P0) ) {

        if ( min(t)>0 )
            warning("analytic solution for proteins is currently untested",
                    "and appears wrong for min(t)>0!",
                    "Use ODE model for correct solutions.")

        ## fix for t>>0 probem: set initial time to 0
        ## TODO: is this correct?
        tp <- t-min(t)
        
        ## total turnover
        gammap <- mu+dp

        ## exponential decay term
        EGP <- exp(-gammap*tp)
    
        ## RNA contributions

        ## beta*int(R(t)):
        ## integrate calculated RNA values
        intr <- rt[,"R"] * exp(gammap * tp)
        ## Trapezoidal rule
        intr <- cumsum((intr[-1] + intr[-nrow(rt)]) / 2) * diff(tp) 

        ## translation rate
        beta <- ell*rho
        intr <- beta * intr
        
        ## protein abundance
        ## TODO: beta *rt[1,"R"] instead of 0
        pt <- EGP*(P0 + c(0, intr))

        ## bind RNA and protein abundances
        rt <- cbind.data.frame(rt, P=pt)
   
    }

    rt
}

## plot labels
## TODO: get axis/unit functions used for stan model script
axis_labels <- c(rmean=expression('\u27E8'*R*'\u27E9'/(n/cell)),
                 ramp=expression(tilde(R)),
                 rampr=expression(tilde(r)),
                 r=expression(R(t)),
                 phi=expression(duty~cycle~phi),
                 tau=expression(period~tau))
