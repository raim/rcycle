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
biplot(gphase, scale=0, xlabs = rep("x", nrow(gphase$x))) 

## WITH SCALING:
plotPC(gphase, arrows=FALSE, scores=TRUE, time.line=TRUE, eigen.axis=TRUE,
       scale=1, arcsinh=FALSE, col=2, pch=19, cex=1)
par(new=TRUE)
bp <- biplot(gphase, scale=1, xlabs = rep("x", nrow(gphase$x)))

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
