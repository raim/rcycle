
library(rcycle) # this library
library(segmenTools) # for num2col only
library(RSpectra) # for eigs_sym

## Compares SVD/PCA on simulated data (sine waves, with noise, just noise)

## 1. generate data genes X samples,
## 2. norm by total counts (single cells),
## 3. row-centering or standardization,
## 4. "Eigengenes": do not transpose
##     - raw data: svd, eigen(cov()), prcomp(scale.=FALSE, center=FALSE)
##     - row-centered: same same,
##     - column centering: svd, eigen(cor), prcom(scale.=FALSE)
##     - column scaling: svd, eigen(cor), prcom(scale.=FALSE)
## 5. single cell: transpose and do SVD=PCA,

## TODO: get sine wave simulation from cohorts_simulate.R and
## PWM simulation from models.R and pwmode.R in ChemostatData/src/
