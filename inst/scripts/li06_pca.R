
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

##rcycle.datafile <- paste0('~/programs/rcycle/data/',did,'_data.rda')
##load(rcycle.datafile)

pdat <- li06_data$data
time <- li06_data$time
tosc <- li06_data$tosc

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
cat(paste('Running eigenvalue analysis:\t',did,' - TAKES LONG\n'))

## eigen(cor(x)) == prcomp(x, scale=TRUE)
gev <- eigen(cor(gdat), symmetric=TRUE) 

## eigen(cor(t(x))) == prcomp(t(x), scale=TRUE); row-norm BEFORE t()
sev <- RSpectra::eigs_sym(cor(sdat), k = ncol(gdat))

## eigen(cor(t(x))) == prcomp(t(x), scale=TRUE); row-norm AFTER t()
## NOTE: cov required, see above scale=FALSE
cev <-  RSpectra::eigs_sym(cov(cdat), k = ncol(gdat)) 

## cov/cor vs. scale
gevv <- eigen(cov(gdat), symmetric=TRUE)
gphasev <- prcomp(gdat, scale.=FALSE)

## column centering vs. scale
gdatc <- t((t(gdat) - apply(gdat,2,mean))/apply(gdat,2,sd))
gphasec <- prcomp(gdatc, scale.=FALSE)

### ANALYZE

## BIPLOT vs. plotPC
## ... should, when multiplied together, approximate X (that's the
## whole point).
## https://stats.stackexchange.com/questions/66926/what-are-the-four-axes-on-pca-biplot
## https://stats.stackexchange.com/questions/141085/positioning-the-arrows-on-a-pca-biplot

## TODO: test scaling more exactly!
plotPC(gphase, arrows=FALSE, scores=TRUE, time.line=TRUE, eigen.axis=TRUE,
       arcsinh=FALSE) # eigenvectors/rotation matrix are samples
par(new=TRUE)
biplot(gphase, scale=0) 

## WITH SCALING:
plotPC(gphase, arrows=FALSE, scores=TRUE, time.line=TRUE, eigen.axis=TRUE,
       scale=1, arcsinh=FALSE, col=2, pch=19, cex=1)
par(new=TRUE)
bp <- biplot(gphase, scale=1)

## scale eigenvectors!
n <- nrow(gphase$x)
lam <- gphase$sdev * sqrt(n)
gphase$scaled <-t(t(gphase$rotation) * lam) 
    
plot(gphase$rotation[,1], gphase$rotation[,2], type='b')
par(new=TRUE)
plot(gphase$scaled[,1], gphase$scaled[,2], type='b', col=2, pch=4, axes=FALSE,
     xlab=NA, ylab=NA)
mtext("gphase$scaled[,1]", 3, par('mgp')[1], col=2)
mtext("gphase$scaled[,2]", 4, par('mgp')[1], col=2)
axis(3, col.axis=2)
axis(4, col.axis=2)
abline(h=0)
abline(v=0)

## minor consequence for phase angle?
plot(atan2(gphase$rotation[,2], gphase$rotation[,1]),
     atan2(gphase$scaled[,2], gphase$scaled[,1]))
abline(a=0, b=1)

## biplots for "conventional" PCA with genes as columns
plotPC(sphase, arrows=FALSE, ccol=2, cpch=19)
lines(sphase$x[,1], sphase$x[,2], col=2) # rotated data are samples
##biplot(sphase) 
plotPC(cphase, arrows=FALSE, ccol=2, cpch=19)
lines(cphase$x[,1], cphase$x[,2], col=2) # rotated data are samples
##biplot(cphase) 

### DIRECT CALCULATION as the eigenvalues/vectors of the cor/cov matrix
## @Schwabe2000: Cov(N^T)=1/(n-1) N * N^T, where N is a row-centered
## geneXcell matrix (m‐by‐n (genes‐by‐cells)), and get the
## eigenvectors W, where Cov(N^T) = W^T * D * W, with the diagonal
## matrix of eigenvalues, where the rows of W^T are the eigenvectors
## of Cov⁡(N^T). The principal components then are P = W^T * N. A row
## vector p_k, of P contains the PC scores (or amplitudes) of all
## cells,

## eigen(cor(x)) == prcomp(x, scale=TRUE)
plot(gev$vectors[,1], gphase$rotation[,1])
abline(a=0, b=1)

## eigen(cor(t(x))) == prcomp(t(x), scale=FALSE); row-norm AFTER t()
plot(sev$vectors[,1], sphase$rotation[,1])
abline(a=0, b=1)

## eigen(cov(t(x))) == prcomp(t(x), scale=FALSE); row-norm BEFORE t()
plot(cev$vectors[,1], cphase$rotation[,1]) #NOTE: negative when using eigen!
abline(a=0, b=1)

## cov/cor vs. scale
plot(gevv$vectors[,1], gphasev$rotation[,1])
abline(a=0,b=1)

## column standardization vs. scale
plot(gphase$rotation[,1], gphasec$rotation[,1])
abline(a=0,b=1)


## relation of PCA sdev to eigen values
## var = sd^2
plot(gev$values[1:length(gphase$sdev)], gphase$sdev^2)
abline(a=0,b=1)
            
plot(sev$values[1:length(sphase$sdev)], sphase$sdev^2)
abline(a=0,b=1)
 
plot(cev$values[1:length(cphase$sdev)], cphase$sdev^2)
abline(a=0,b=1)

## COMPARE PSEUDOPHASES
## pseudophases are very similar between approaches
plot(gtheta, stheta)
plot(gtheta, ctheta)
plot(stheta, ctheta)


## EQUIVALENCE BETWEEN eigen and eigs_sym?
## NOTE: while numerical differences exist,
## only the final eigenvalue seems to differ,
## but some errors are present
test.rspectra <- FALSE # SLOW
if ( test.rspectra ) {

    sev2 <- eigen(cor(sdat), symmetric=TRUE)
    cev2 <- eigen(cov(cdat), symmetric=TRUE)

    cors <- sapply(1:ncol(gdat), function(x)
        cor(cev2$vectors[,x], cev$vectors[,x]))

    barplot(log10(abs(cors)))

    errs.idx <- which(abs(cors) < 1)
    plot(c(cev2$vectors[,1:ncol(gdat)]), c(cev$vectors), col=NA)
    for ( i in errs.idx )
        lines(cev2$vectors[,i], cev$vectors[,i], col=i)
    legend('top', legend=errs.idx, lty=1,
           col=errs.idx)
}

### IMPORTANT POINTS:
## * a correct biplot is defined such that multiplication yields
##   would yield the original data!
## * cov/cor vs. scale
## eigen(cov(x)) == prcomp(x, scale=FALSE), and
## eigen(cor(x)) == prcomp(x, scale=TRUE);
## cor is the cov(x,y)/(var(x)*var(y))
## * scale=TRUE is equivalent to column centering!

## * The current pseudophase implementation in rcycle does
##   it unconventionally, where we (a) both row-normalize and standardize
##   columns, and (b) use CELLS as variables (columns) and get pseudophase
##   from the eigenvectors 1 vs. 2, instead of the more conventional approach.
