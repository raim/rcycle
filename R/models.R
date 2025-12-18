
## MODELS for PERIODIC GENE EXPRESSION

## PULSE WAVE GENERATING FUNCTIONS

pwm_simple <-  function(time, tau, phi, theta = 0) {
  # Returns integer 0/1 vector: 1 when pulse is on.
  # theta in radians; tau period in same time units as time.
  omega <- 2 * pi / tau
  # convert phase shift theta (radians) to time shift
  t_shifted <- (time - theta / omega) %% tau
  as.integer(t_shifted < (phi * tau))
}


## pulse wave function - direct and very SLOW calculation
## wikipedia: "Note that, for symmetry, the starting time (t = 0)
## in this expansion is halfway through the first pulse."
## TODO: why do we need to scale?
#' alpha controls the smoothing by exponentially damping higher harmonics; alpha=0: no smoothing; alpha<=0.5: mild smoothing; alpha>0.5: noticable smoothing; alpha>2: towards sinus
#' 
#' t <- seq(0,10, length=100)
#' plot(t, pwf(t=t, tau=1.5, k=1, phi=.4), type='l', col=2)
#' lines(t, pwf(t=t, tau=1.5, k=1, phi=.6, theta=pi), type='l', col=4)
#'
pw_fourier <- function(t=0, k, phi, tau, N=1e4, shift=0, theta=0,
                       start.on=FALSE, alpha=0) {

    ## shift to start with on phase
    ## NOTE: to start at on we'd shift in the cos function:
    ## cos(n*omega*t - n*pi*phi)
    if ( start.on )
        t <- t - phi*tau/2
 
    
    ## shift time to get phase shift (passed in unit of time!)
    t <- t - shift

    ## sum expression
    ## TODO: faster vectorization instead of loop here but we have
    ## both n and t vectors?
    omega <- 2*pi/tau
    sincos <- 0
    for ( n in 1:N )  
        sincos <- sincos +
            1/n * sin(pi*n*phi) * cos(n*omega*t - n*theta)* exp(-alpha*n)
    x <- 2/pi * sincos
    
    ## shift to 0 and k
    x <- k*(phi + x) 
    x
}

## pulse wave via sawtooth: subtract phase-shifted
## NOTE: sawtooth starts with off phase, but we can shift time
## here to ensure it starts with an on phase
## TODO: different formulation to avoid min/max scaling?
pw_sawtooth <- function(t=0, k=1, tau, thoc=.5, alpha=0, shift=0,
                        start.on=FALSE, start.like.pw=TRUE, scale=TRUE) {

    if ( start.like.pw ) # shift to start like pwf with half pulse!
        t <- t - thoc/2
    else if ( start.on ) # shift to start with on phase
        t <- t + (tau-thoc)
        
    ## shift time to get phase shift (passed in unit of time!)
    t <- t - shift
    
    st <-  atan(pracma::cot(pi* t      /tau)) 
    st2 <- atan(pracma::cot(pi*(t+thoc)/tau))

    ## TODO: implement smoothing of saw tooth function!
    if ( alpha!=0 ) {
        stop('smoothing not implemented yet for saw tooth function')
        sig2 <- (t+thoc)/alpha
        sig2 <- 1/(1-exp(-sig2))
        sig1 <- t/alpha
        sig1 <- 1/(1-exp(-sig1))
        
        st2 <- st2*sig2
        st <- st*sig1
        ##return(sig)
    } 
    
    x <- 1/pi * (st2 - st)
    
    ## scale between 0 and k
    ## TODO: use different formula to avoid min/max scaling?
    if ( scale )
        x <- k*(x-min(x))/(max(x)-min(x)) # x <- k*thoc/tau + x
    x
}

### ODE Models


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
get_rmean <- function(k, gamma, k0, dr, mu, phi, tau,
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

    unname(rmn)
    
}


#' Calculate abundance amplitudes from rates and times.
#'@export
get_ramp <- function(gamma, dr, mu, phi, tau, relative = TRUE, k, k0,
                     model = c('k', 'dr', 'k_dr', 'k_dr_k0'), ...) {

    if ( model %in% c('dr', 'k_dr', 'k_dr_k0') ) {

        ## calculate relative amplitude
        if ( missing(gamma) )
            gamma <- dr+mu

        gt <- gamma*tau

        ## TODO: use expm1 version?
        term2 <- phi/2 * pracma::coth(gt*(1-phi)/2)

        term1 <- switch (model,
                         k_dr = 1/gt,
                         dr = (1-1/phi)/gt,
                         k_dr_k0 = (1-k0/(k*phi))/gt)

        ramp <- 1/(term1 + term2)

        if ( missing(k0) ) k0 <- NA

        if ( !relative ) {
            rmean <- get_rmean(
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

    unname(ramp)
    
}

get_rna <- function() {

    ## TODO: wrapper for get_rmean and get_ramp
}



get_times <- function(model = c('k', 'dr', 'k_dr', 'k_dr_k0'),
                      k, dr, k0, lower = 1e-6, upper = 20, ...) {
}

### INTERNAL ROOT FINDING FUNCTIONS



get_rates <- function(model = c('k', 'dr', 'k_dr', 'k_dr_k0'),
                      a = NA, A = NA, R = NA, Rmin = NA, Rmax =NA,
                      phi = NA, tau = NA, mu = NA,
                      k = NA, k0 = NA, gamma = NA, 
                      lower = 1e-6, upper = 20, tol = 1e-9,
                      verb = 0, ...) {

    ## relative amplitude is required for all
    if ( all(is.na(a)) ) a <- A/R

    if ( !model %in% c('k')  ) {
        
        ## get transcription rate from abs. amplitude
        if ( all(is.na(k)) ) {
            if ( all(is.na(A)) )
                warning('transcription rate k requires absolute amplitude A')
            k <- get_transcription(A = A, phi = phi, tau = tau, model = model)
        }
        ## only required for k_dr_k0
        ## get Rmin from Rmax and k
        if ( is.na(Rmin) & !is.na(Rmax) ) Rmin <- Rmax - k*phi*tau
    }
    
    dr <- unlist(Map(get_degradation,
                     model = model,
                     a=a, A=A, R=R, Rmin=Rmin, Rmax=Rmax,
                     phi = phi, tau = tau, mu = mu,
                     k = k, 
                     lower = lower, upper = upper, tol = tol,
                     verb = verb)   )
    
    if ( model %in% c('k') ) {

        gamma <- dr
        if ( !is.na(mu) ) gamma <- dr + mu
   
        k <- get_transcription(R = R, gamma = gamma, phi = phi, model = model)

        if ( length(k)==1 )
            k <- rep(k, length(dr))
    }

    ## add k0
    
    res <- cbind(k=k, dr=dr)
    rownames(res) <- NULL
    res
}

#' Calculate transcription rate from average abundance and degradation rates.
#' from base equation for average abundance
#' TODO: map fluorescence to a rough transcript/cell number
#' @param R average abundance
#' @param A absolute amplitude, for biphasic models
get_transcription <- function(R = NA, A = NA, phi = NA,
                              dr = NA, mu = NA, gamma = NA, tau = NA,
                              model = c('k', 'dr','k_dr', 'k_dr_k0')) {

    if ( model %in% c('k') ) {
        if ( !is.na(dr) & !is.na(mu) )
            gamma <- mu+dr
        k <- R*gamma/phi
    } else 
        k <- A/(phi*tau)
    
    k
}

get_basal <- function(k, gamma, tau, phi, Rmin, Rmax, verb = 1) {

    ##beta <- exp(gamma*tau*(1-phi))
    betam1 <- expm1(gamma*tau*(1-phi)) 
    
    if ( !missing(Rmin) ) 
        k0min <- k0 <- (Rmin - k*phi*tau/betam1)*gamma
    if ( !missing(Rmax) ) 
        k0max <- k0 <- (Rmax - k*phi*tau*(1 + 1/betam1))*gamma
    if ( exists('k0min', mode = 'numeric') &
         exists('k0max', mode = 'numeric') ) {
        if ( verb>0 ) cat(paste('mean of rmin- and rmax-based basal rates\n'))
        k0 <- (k0min+k0max)/2    
    }
    k0
}


#' Fit a transcript degradation rate from relative amplitude and
#' oscillation parameters, using the \code{\link[stats]{uniroot}}  function
#'
#' @param a relative abundance amplitude (max(x)-min(x)/mean(x)).
#' @param R mean abundance, required for model with basal transcription.
#' @param Rmin minimal abundance, required for model with basal transcription.
#' @param k transcription rate, required for model with basal transcription.
#' @param phi duty cycle.
#' @param tau period.
#' @param mu growth rate, used only for correction (gamma=growth+degradation).
#' @param lower the lower end point of the interval to be  searched by \code{\link[stats]{uniroot}}.
#' @param upper the upper end point of the interval to be  searched by \code{\link[stats]{uniroot}}.
#' @param model 
#' @param tol error tolerance of the \code{\link[stats]{uniroot}}  call.
#' @param ... further parameters to \code{\link[stats]{uniroot}} 
get_degradation <- function(a,  R, Rmin, Rmax, k,
                            phi, tau, mu,
                            model, ## must exist as root finding function
                            lower = 1e-6, upper = 1e4, tol = 1e-9,
                            verb = 0, ...) {

    

    
    ## get model-specific function for root-finding,
    ## each returns x = gamma * tau
    rootf <- get(paste0('root_', model), mode = 'function')
    

    ## Solve f(x) = 0 for x in a reasonable range
    solution <- try(stats::uniroot(rootf, a=a, phi=phi, R=R, Rmin=Rmin,
                                   lower = lower, upper = upper, tol = tol),
                    silent = verb==0)
    if ( class(solution)=="try-error" ) {
        if ( verb>0 )
            cat(paste0('Model <', model, '> failed with:',
                       '\n\ta=', a,
                       '\n\tR=', R,
                       '\n\tRmin=', Rmin,
                       '\n\tRmax=', Rmax,
                       '\n\tk=', k,
                       '\n\tphi=', phi,
                       '\n\ttau=', tau,
                       '\n\tmu=',  mu, '\n'))
        return(NA)
    }
        
    ## Recover gamma from x = gamma * tau
    x_sol <- solution$root
    gamma <- x_sol / tau

    ## correct for growth
    ## gamma = degradation + growth
    if ( !missing(mu) )
        gamma <- gamma - mu

    gamma
}

## ROOT FINDING FUNCTIONS USED IN get_degradation

## transcription with const. degradation,
## where x = gamma * tau
root_k <- function(x, a, phi, R = NA, Rmin = NA) {
    lhs <- a * phi
    ##term1 <- (1 - exp(-phi * x)) / (1 - exp(-x))
    ##term2 <- (1 - exp(x * (phi - 1)))
    term1 <- (- expm1(-phi * x)) / (- expm1(-x))
    term2 <- ( - expm1(x * (phi - 1)))
    rhs <- term1 * term2
    return(lhs - rhs)
}
    


## anti-phasic transcription/degradation,
## where x = gamma * tau
root_k_dr_coth <- function(x, a, phi, R = NA, Rmin = NA) {
    
    lhs <- 1/a
    term1 <- 1/x
    term2 <- phi/2 * pracma::coth((x-x*phi)/2)
    rhs <- term1+term2
    return(lhs - rhs)
}
root_k_dr <- function(x, a, phi, R = NA, Rmin = NA) {
        
    lhs <- 1/a
    term1 <- phi/expm1(x * (1-phi)) # = exp() + 1, where 1 cancels below
    term2 <- phi/2
    term3 <- 1/x
    rhs <- term1 + term2 + term3
    return(lhs - rhs)
}

## degradation pulse
## where x = gamma * tau
root_dr <- function(x, a, phi, R = NA, Rmin = NA) {
    
    lhs <- 1/a
    term1 <- (1-1/phi)/x
    term2 <- phi/2 * pracma::coth((x-x*phi)/2)
    rhs <- term1+term2

    return(lhs - rhs)
}


## anti-phasic transcription/degradation and basal transcription
## where x = gamma * tau
root_k_dr_k0_coth <- function(x, a, phi, R, Rmin) {

    A <- a*R
    lhs <- (R - Rmin)/A

    term1 <- 1/x
    term2 <- 1/expm1(x*(phi-1))
    term3 <- phi/2 * pracma::coth((x-x*phi)/2)
    rhs <- term1 - term2 + term3
    return(lhs - rhs)
}
root_k_dr_k0 <- function(x, a, phi, R, Rmin) {
    
    lhs <- R - Rmin
    betam1 <- expm1(x*(1-phi)) # exp() +1
    term1 <- (phi-1) /betam1 
    term2 <- phi/2
    term3 <- 1/x 
    rhs <- a*R*(term1 + term2 + term3)
    ##rhs <- phi*k*tau*(term1 + term2 + term3)
    return(lhs - rhs)
}


### ANALYTIC MODEL

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
#' @param tau Oscillation period.
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
                  phi, tau,
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
    omega <- 2*pi/tau

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

### SOME PLOT UTILS

## plot labels
## TODO: get axis/unit functions used for stan model script
axis_labels <- c(rmean=expression('\u27E8'*R*'\u27E9'/(n/cell)),
                 ramp=expression(tilde(R)),
                 rampr=expression(tilde(r)),
                 r=expression(R(t)/(n/cell)),
                 phi=expression(duty~cycle~phi),
                 tau=expression(period~tau/h),
                 k=expression(k/(n/h)),
                 dr=expression(delta[R]/h^-1),
                 gamma=expression(gamma/h^-1))
