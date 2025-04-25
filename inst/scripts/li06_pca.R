
### understanding PCA for raim, 20250411

## * TESTING PCA with genes or samples as variables,
## * comparing plotPC with biplot, where values are scaled
##   bei eigenvalues such that the product would yield the original matrix.
## * trying to align with methods part of @Schwabe2000, and understand
##   how eigen(cor/cov) and prcomp(scale=TRUE/FALSE) are related,
##   and what role column and row normalization play.

## * TODO:
## * test the same for simulated, single cell data, and cohort approach!
## * convert this to an Rmarkdown file, perhaps as a vignette for rcycle!
## * compare e.g. to Monocle and diffusion time approaches, and to WGCNA. 

library(rcycle) # this library
library(segmenTools) # for num2col only
library(RSpectra) # for eigs_sym

## use latest during development
source("~/programs/rcycle/R/phases.R")
source("~/programs/rcycle/R/plot.R")

did <- 'li06'
data('li06_data')
##ddid <- paste0(did,'_data')
##li06_data <- get(ddid)

rcycle.datafile <- paste0('~/programs/rcycle/data/',did,'_data.rda')
load(rcycle.datafile)

pdat <- li06_data$data
time <- li06_data$time
tosc <- li06_data$tosc

states <- li06_data$states
scol <- li06_data$col

## TODO: separate script, comparing
## genes and states
use.states <- FALSE
if ( use.states ) pdat <- states

colnames(pdat) <- sub('li06_', '', colnames(pdat))

## time point coloring by oscillation cycle
otime <- time - min(time)
while(any(otime>tosc))
    otime[otime>tosc] <- otime[otime>tosc] - tosc
tcol <- segmenTools::num2col(otime)


### GENE-CENTERED
## VARIABLES/COLUMNS are samples/cells, OBSERVATIONS/ROWS are genes
## row-normalized and UNTRANSPOSED
cat(paste('Running PCA on genes:\t',did,'\n'))
gdat <- pdat - apply(pdat, 1, mean)
gphase <- prcomp(gdat, scale.=TRUE)

### SAMPLE-CENTERED with genes as columns!
## VARIABLES/COLUMNS are genes, OBSERVATIONS/ROWS are samples/cells

## row-norm over cells
cat(paste('Running PCA on samples:\t',did,'with prior row-norm over samples\n'))
sdat <- t(pdat) # transpose BEFORE row-norm 
sdat <- sdat - apply(sdat, 1, mean)  # row-norm over cells!
sphase <- prcomp(sdat, scale.=TRUE)

## row-norm over genes: as in @Schwabe2000?
cat(paste('Running PCA on samples:\t',did,'with prior row-norm over genes\n'))
cdat <- pdat - apply(pdat, 1, mean)  # row-norm over genes
cdat <- t(cdat) # transpose AFTER row-norm 
cphase <- prcomp(cdat, scale.=FALSE)  # NOTE: -> colnorm!

## pseudophases for all approaches
cat(paste('getting pseudophases:\t',did,'\n'))
gtheta <- atan2(gphase$rotation[,2], gphase$rotation[,1])
stheta <- atan2(sphase$x[,2], sphase$x[,1])
ctheta <- atan2(cphase$x[,2], cphase$x[,1])


## Calculate eigenvalues of the cor/cov matrices directly
cat(paste('Running eigenvalue analysis:\t',did,' - TAKES LONGER\n'))

## eigen(cor(x)) == prcomp(x, scale=TRUE)
gev <- eigen(cor(gdat), symmetric=TRUE) 

## eigen(cor(t(x))) == prcomp(t(x), scale=TRUE); row-norm BEFORE t()
sev <- RSpectra::eigs_sym(cor(sdat), k = ncol(gphase$x))

## eigen(cov(t(x))) == prcomp(t(x), scale=FALSE); 
## NOTE: cov required, see above scale=FALSE
cev <-  RSpectra::eigs_sym(cov(cdat), k = ncol(gphase$x)) 

## cov/cor vs. scale
gevv <- eigen(cov(gdat), symmetric=TRUE)
gphasev <- prcomp(gdat, scale.=FALSE)

## column centering vs. scale
gdatc <- t((t(gdat) - apply(gdat,2,mean))/apply(gdat,2,sd))
gphasec <- prcomp(gdatc, scale.=FALSE)

## TODO: prcomp vs. svd vs. eigen and ALL vs. row-centering!
## 1. svd, eigen(cov), prcomp(scale.=FALSE, center=FALSE),
## 2. additional row-centering

psv <- svd(pdat) # raw
gsv <- svd(gdat) # row-norm over genes
ssv <- svd(sdat) # transposed and row-norm over cells
csv <- svd(cdat) # row-norm and transposed (i.e. col-norm)



### TODO: COMPARE with pca of cohort states
