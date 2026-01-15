
library(rcycle)

##. temp requirements
library(viridis)
library(segmenTools)
source('/home/raim/programs/rcycle/R/models.R')


out.path <- '/home/raim/programs/rcycle/vignettes'
ftyp <- 'pdf'

## TODO:
## * load fitted parameters for RP and histone genes,
## * solve phi(mu) for ribosomes and histones,
## * constant duration: histones: S phase, ribosomes: budding?



## parameters for average RP gene 

## transcription
dr <- 2.799240  # RNA degradation rate, 1/h
k <- 2.638813e+02 # transcription rate, n/h

# translation
ell <- 250 # translation rate, n/h
rho <- 3 # ribosomes/RNA
dp <- 6.953300e-02 # protein degradation rate, 1/h

## growth rates 
mu <- seq(0,.5,length.out=100)

### RNA vs. PROTEIN ABUNDANCE vs. GROWTH RATE

mus <- seq(0,.35,.01)
mlim <- c(0,.4)
phis <- rev(seq(0,1,.05))

## phi colors: pulse wave, smoothed and half-width and period switching
kcols <- c(viridis::viridis(3)[c(1,3,2)],"#00B7EBFF")
kcols[2] <- "#FF4D01FF"
phi.cols <- colorRampPalette(kcols[c(2,3)])(length(phis)) 

## R/P for different growth rates and phihoc
rmat <- pmat <- rp <- matrix(NA, nrow=length(phis), ncol=length(mus))
for ( i in 1:length(phis) ) {
    rmat[i,] <- get_rmean(dr = dr, k = k, phi = phis[i], mu = mus, model = 'k')
    pmat[i,] <- get_pmean(R = rmat[i,], rho = rho, l = ell, dp = dp, mu = mus)
    rp[i,] <- pmat[i,]/rmat[i,]
}

## 60000 nucleosomes per genome (12e6/150), 2 histones per nucleosome,
## 2 genomes per cell
P_target <- 2*2*6e4 
## average nucleosome depends on fraction of budding cells
tau_bud <- rep(1, length(mus))
phi_bud <- tau_bud*mus/log(2)
## mean histone concentration in population
P_pop <- P_target*(1+phi_bud)

for ( i in 1:4 ) {

    pcols <- phi.cols
    if ( i==1 ) {
        pcols[] <- NA
        pcols[length(pcols)/2] <- phi.cols[length(pcols)/2]
    }
    
    plotdev(file.path(out.path, paste0("growthrate_rna_", i)),
            type=ftyp, width=2.5, height=2.5, res=300)
    par(mfrow=c(1,1), mai=c(.5,.5,.1,.1),
        mgp=c(1.3,.3,0), tcl=-.25)
    matplot(mus, t(rmat), type="l", lty=1,
            col=pcols, lwd=2, xlim=mlim,
            xlab=expression("growth rate"~mu/h^-1),
        ylab=expression(transcripts/(n/cell)))

    ## indicate phi
    if ( i>1 ) {
        arrows(x0=.37, y0=0, y1=85, lwd=3, length=.1, col=pcols[1])
        nl <- length(pcols)-1; ny <- 85-.07
        shadowtext(rep(.36,2), c(.1,84), c("0","1"), pos=4,
                   col=pcols[c(length(phis),1)])
        for ( j in nl:1 )
            lines(x=rep(.37,2), y=.07+c(j*ny/nl, (j-1)*ny/nl),
                  col=rev(pcols)[j], lwd=4)
        shadowtext(x=.35, y=40, expression(varphi),
                   pos=4, col=1, cex=1.2)
    }
    dev.off()
   
    plotdev(file.path(out.path, paste0("growthrate_protein_",i)),
            type=ftyp, width=2.5, height=2.5, res=300)
    par(mfrow=c(1,1), mai=c(.5,.5,.1,.1),
        mgp=c(1.3,.3,0), tcl=-.25, yaxs="i")
    matplot(mus, t(pmat)/1e5, type="l", lty=1, ylim=c(0,3.5),
            col=pcols, lwd=1.5, xlim=mlim,
            xlab=expression("growth rate"~mu/h^-1),
            ylab=NA)
    mtext(expression(proteins/(10^5%*%n/cell)), 2, 1.1)

    ## indicate example proteins
    if ( i>2 ) {
        ##legend("topright", " - -   protein\nhomeostasis", pch=NA, 
        ##       box.col=NA, bg="#ffffff00", seg.len=0, lty=NA, lwd=1,
        ##       text.col=2, text.font=2)

        ## ribosomes
        arrows(x0=0.03, x1=0.3, y0=1, y1=2, length=.1,
               lwd=5, col="#ffffffbb", code=0)
        arrows(x0=0.03, x1=0.3, y0=1, y1=2, length=.1,
               lwd=2, col=2, code=0, lty=1)
        ## TF
        lines(mus, pmat[nrow(pmat)-3,]/1e5, type="l", col=5, lwd=2)

        ## histones: once wrong, and once more correctly!
        if ( i == 4 ) {
            arrows(x0=0.03, x1=0.33, y0=.6, y1=.6, length=.1,
                   lwd=5, col="#ffffffbb", code=0)
            arrows(x0=0.03, x1=0.33, y0=.6, y1=.6, length=.1,
                   lwd=2, col=1, code=0, lty=1)
            hlab <- 'histones?'
        } else {
            hlab <- 'histones!'
            lines(mus, P_pop/1e5, col = "#ffffffbb", lwd=5)
            lines(mus, P_pop/1e5, col = 1, lwd=2)
        }
        legend("topright",c("ribosomes","e.g. LacI",hlab),
               col=c(2,5,1),
               seg.len=.75, lty=1, lwd=2, bg="#ffffff77", box.col=NA)

    }
    if ( i > 2 ) {
        arrows(x0=.37, y0=.07, y1=1.5, lwd=3, length=.1, col=pcols[1])
        nl <- length(pcols)-1; ny <- 1.5-.07
        shadowtext(rep(.36,2), c(.1,1.45), c("0","1"), pos=4,
                   col=pcols[c(length(phis),1)])
        for ( j in nl:1 )
            lines(x=rep(.37,2), y=.07+c(i*ny/nl, (j-1)*ny/nl),
                  col=rev(pcols)[j], lwd=4)
        shadowtext(x=.35, y=.75, expression(varphi),
                   pos=4, col=1, cex=1.2)
    }
    dev.off()
}

